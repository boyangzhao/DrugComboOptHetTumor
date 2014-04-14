%% Simulation Class for combination therapy optimization on heterogeneous populations
% Author: Boyang Zhao
% Initiates simulation, generates population structures, and calls optimization

classdef (Sealed) Simulation < handle & OptimizationModels & Analysis
    
    %% Public methods
    methods (Access = public)
        function m = Simulation(dsets)
            m = m@OptimizationModels(dsets);
            m = m@Analysis(dsets);
        end
        
% SIMULATION
        function simulate(m)
            %Determine optimal therapies for given/simulated population structures
            
            %check for setting errors
            if(m.drugsN > 0 && length(m.drugsInclude) > m.drugsN)
                error('SettingsError:DrugsConstraint', ...
                      'The constraint on the optimal number of drugs is smaller than the numbers of drugs to include');
            end
            
            %check optimization model to use
            if(isempty(m.optModel))
                disp('No optimization model specified');
            end
            
            m.optModelName = m.optModelMap(m.optModel);
            disp(['Dataset: ',m.dName]);
            disp(['Optimization model: ',m.optModelName]);
            
            %check therapy settings
            if(m.sametox)
                disp('Using same toxicity profile for all drugs');
                m.dTox = ones(size(m.dTox)); %change toxicity to all ones
            end
            
            %start simulation initialization
            %enumerate population structures (i.e. subpopulations and their proportions)
            popstructs = m.enumeratePopStruct(); %popstruct is abbreviated version
            
            disp('**********************');
            disp('Starting simulation...');
            
            %Note: can preallocate m.results for performance
            %optimize for each popstruct
            for p = 1:length(popstructs)
                %derive the expanded popstruct (sparse vector)
                m.popstruct_e = m.expandPopstruct(popstructs(p).subpops, popstructs(p).contris);
                
                %static optimization
                m.optimize(); %results stored in m.solutions
                
                %update progress
                if Utilities.comprval(rem(p,length(popstructs)/10),0)
                    disp([num2str(p*100/length(popstructs),'%6.0f') '% Completed']);
                end
                
                 m.results = [m.results;m.solutions];
            end
        end
    end
    
    %% Private methods
    methods (Access = private)    
% GENERATE POPULATION STRUCTURES
        function popstructs = enumeratePopStruct(m)
            %Generate a list of population structures (abbreviated version)
            %each popstruct has the following fields: N, subpops, contris
            %If m.popStructGiven is not empty, will use popstruct from here instead
            %calls either generateEnumeratedPopStructs() or generateRandPopStruct()
            
            if(m.popStructGiven(1).N ~= 0)
                disp('Using given population structure...');
                %will ignore m.montecarlo, m.simulateN, m.popStructSet, m.popStructSetPrior, and m.subpopN
                %unless contris does not add up to 1, then the remaining proportions will be Monte Carlo sampled if 
                %m.monetecarlo is set to true; otherwise will return error;
                %this will be checked 
                popstructs = m.popStructGiven;
                
                if(m.popFillContris)
                    %will check through all the popstructs, and fill sample and fill in the rest of the popstruct 
                    %if contris does not add up to use (for this m.montecarlo must be set to true)
                    disp('m.popFillContris is set to true; will Monte Carlo sample to fill the rest if needed...');
                    if(~m.montecarlo | m.simulateN < 1)
                        error('SettingsError:MonteCarlo', ...
                            'If m.popFillContris is set to true, m.montecarlo need to be set to true and relevant settings defined (m.popStructSet and m.popStructSetPrior)');
                    end
                    
                    popsize = length(popstructs);
                    popfilled = 0;
                    
                    %save current popStructSet; m.popStructSet will be temporarily 
                    %updated (if filling is needed) in loop below
                    popStructSet_original = m.popStructSet;
                    for i = 1:popsize                        
                        if(sum(popstructs(i).contris) < 1)
                            %update popStructSetPrior and popStructSet
                            %with the current fixed subpopulations excluded
                            [~,selectfixed,~] =setxor(popStructSet_original, popstructs(i).subpops);
                            m.popStructSet = popStructSet_original(selectfixed);
                            m.popStructSetPrior = (1/length(m.popStructSet))*ones(size(m.popStructSet));
                            contrisTotalLeft = 1-sum(popstructs(i).contris);
                            workingpop = popstructs(i);
                            for j = 1:m.simulateN
                                %draw a random mixed set of subpopulations for the remainder of the population
                                p_remainder = m.generateRandPopStruct();
                                %drop the ones that are less than the minimum proportion (set to 1%)
                                mincontri = 1;
                                selectsubpop = ~( (contrisTotalLeft.*p_remainder.contris) < (mincontri/100) );
                                p_remainder.N = sum(selectsubpop);
                                p_remainder.subpops = p_remainder.subpops(selectsubpop);
                                p_remainder.contris = p_remainder.contris(selectsubpop);
                                
                                %combine the fixed subpopulations and the remaining random set of subpopulations together
                                p = workingpop;
                                p.N = p.N + p_remainder.N;
                                p.subpops = [p.subpops p_remainder.subpops];
                                p.contris = [p.contris contrisTotalLeft.*p_remainder.contris];
                                if(j == 1)
                                    popstructs(i) = p; %replace the popstruct with the popstruct now with remainder subpops filled in
                                else
                                    popstructs = [popstructs;p]; %for any additional popstruct, add it to end of list
                                end
                            end
                            
                            popfilled = popfilled + 1;
                        end
                    end
                    
                    m.popStructSet = popStructSet_original;
                    
                    disp(['Number of populations given: ' num2str(popsize)]);
                    disp(['Number of populations require to be filled: ' num2str(popfilled)]);
                    disp(['For each population to be filled, number of populations randomly drawn: ' num2str(m.simulateN)]);
                end
                
            else
                %Generate list of popstructs with either Monte Carlo or enumeration
                
                %initialize popstructs
                popstructs = m.getEmptyPopstruct(m.simulateN);
                
                if(m.montecarlo)
                    disp('Starting Monte Carlo generation of population structures...');
                    if(isempty(m.popStructSetPrior))
                        disp('Probabilities for choosing subpopulations is not given, defaulting to uniform distribution...');
                        m.popStructSetPrior = (1/length(m.popStructSet))*ones(size(m.popStructSet));
                    else
                        disp('Using the given probabilities for choosing subpopulations...');
                    end
                    
                    for i = 1:m.simulateN
                        popstructs(i) = m.generateRandPopStruct();
                    end
                else
                    disp('Starting population enumeration...');
                    popstructs = m.generateEnumeratedPopStructs();
                end
            end
        end
        
        function popstructs = generateEnumeratedPopStructs(m)
            %Generate an eumerated list of population structures (abbreviated version)
            %will only generate for m.subpopN is equal to 2
            
            mincontri = 1; %minimum subpopulation proportion, in %
            popstructs = [];
            if(~Utilities.comprval(m.subpopN,-1))
                
                if(Utilities.comprval(m.subpopN,2))
                    popstruct = struct('N',[],'subpops',[],'contris',[]);
                    
                    combos = combnk(1:m.subpopTotal, m.subpopN);
                    a = (mincontri:100-mincontri)';
                    subpopcontris = [a flipud(a)];
                    subpopcontris = subpopcontris/100;
                    
                    popstructs = repmat(popstruct,size(combos,1)*size(subpopcontris,1),1);
                    for i = 1:length(combos)
                        for j = 1:size(subpopcontris,1)
                            popstructs((i-1)*size(subpopcontris,1)+j) =  struct('N',m.subpopN,'subpops',combos(i,:),'contris',subpopcontris(j,:));
                        end
                    end
                else
                    disp('The population size has to be 2 for the eumeration option.');
                end
            end
        end
        
        function popstruct = generateRandPopStruct(m,mincontri)
            %Generate a random population structure (abbreviated version)
            %popstruct has the following fields: N, subpops, contris ([0 1])
            %REQUIRES: m.popStructSet
            %                   m.popStructSetPrior
            %                   m.subpopN [optional]
            
            if nargin < 2
                mincontri = 1; %minimum subpopulation proportion, in %
            end
            
            %pick number of subpopulations
            if(Utilities.comprval(m.subpopN,-1))
                subpopN = randi(length(m.popStructSet),1);
            else
                subpopN = m.subpopN;
            end
            
            %pick subpopulations without replacement
            pset = m.popStructSet;
            pset_prior = m.popStructSetPrior;
            subpops = zeros(1,subpopN);
            
            for n = 1:subpopN
                if length(pset) > 1
                    subpop = randsample(pset,1,true,pset_prior);
                else
                    subpop = pset;
                end
                subpops(n) = subpop;
                pset_prior = pset_prior(pset ~= subpop);
                pset = pset(pset ~= subpop);
            end
            
            %pick contributions (in unit of % first; then convert this to decimal)
            subpopcontris = ones(1,subpopN)*mincontri;
            v_substract = subpopN*mincontri;
            for i=1:subpopN-1
                v = randi([0 100-v_substract]);
                v_substract = v_substract+v;
                subpopcontris(i) = subpopcontris(i)+v;
            end
            subpopcontris(subpopN) = subpopcontris(subpopN)+100-v_substract;
            subpopcontris = subpopcontris/100.0;
            
            popstruct = struct('N',subpopN,'subpops',subpops,'contris',subpopcontris);
        end
        
    end
    
    %% Simulation specific utility methods
    methods (Access = public)
        function pop_e = expandPopstruct(m, subpops, contris)
            %Creates an expanded sparse form of popstruct (abbreviated version) specified by subpopulations and their contributions
            %OUTPUT: pop_e: final output is a m.subpopTotal by 1 vector; length of pop_e = total number of subpopulations
            pop_e = zeros(1,m.subpopTotal);
            pop_e(subpops) = contris;
            pop_e = pop_e';
        end
        
        function [subpops contris] = abbrPopstruct(m, pop_e)
            %Creates an abbreviated form of popstruct
            %INPUT: pop_e: final output is a 1 by m.subpopTotal vector
            %OUTPUT: subpops: final output is a 1 by m.subpopTotal vector
            %               contris: corresponding subpopulation proportions ([0 1])
            pop_e = pop_e';
            [~, subpops, contris] = find(pop_e);
        end
        
        function p = readinPopstruct(m, resultsFilename)
            %Reads in popstruct from specified results file and output a list of popstructs
            %INPUT: simulation results struct (popstruct is in expanded sparse version)
            %OUTPUT: popstruct in abbreviated version (size: n by 1 structs)
            disp('Reading popstruct from file...');
            r = load(resultsFilename,'m');
            resultsN = size(r.m.results,1);
            p = m.getEmptyPopstruct(resultsN);
            
            for i = 1:resultsN
                [subpops, contris]  = m.abbrPopstruct(r.m.results(i).pop);
                N = length(subpops);
                p(i) = struct('N',N,'subpops',subpops,'contris',contris);
            end
        end
        
        function p = getEmptyPopstruct(m, n)
            %Return an array of empty popstructs, n by 1
            %default n = 1
            if(nargin <= 1), n = 1; end;
            
            ep = struct('N',0,'subpops',[],'contris',[]);
            p = repmat(ep,n,1);
        end
    end
    
end


%% Modeling Class for combination therapy optimization on heterogeneous populations
% Author: Boyang Zhao
% Version: 2.10

classdef Modeling < handle
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Input database
    properties (SetAccess = private, GetAccess = public)
        dName;
        drug_labels, subpop_labels;
        dEff;
        dTox_labels;
    end
    
    properties (SetAccess = protected, GetAccess = public)
        dTox;
    end
    
    %Input settings
    properties (SetAccess = public, GetAccess = public)
        outputsdir;
        subpopN, drugsN, sametox, dToxMax, drugsInclude, changeToSubpopOnly;
        montecarlo, simulateN, popStructSet, popStructSetPrior, popFillContris;
        listAllSolutions;
        popStructGiven, therapyToCompare;
    end
    
    %Derived parameter (based on input parameters)
    properties (SetAccess = private, GetAccess = public)
        drugsTotal, subpopTotal;
    end
    
    %Internal variables
    properties (SetAccess = protected, GetAccess = protected)
        popstruct_e;
    end
    
    properties (SetAccess = protected, GetAccess = public)
        results; %vector of m.solutions
        resultsCombined;
    end

    methods (Access = public)
        function m = Modeling(dsets)
            %given inputs
            m.dName = dsets.dName;
            m.drug_labels = dsets.drug_labels;
            m.subpop_labels = dsets.subpop_labels;
            m.dEff = dsets.dEff;
            m.dTox = dsets.dTox;
            m.dTox_labels = dsets.dTox_labels;
            
            %default settings
            m.montecarlo = true;
            m.simulateN = 1000;
            m.subpopN = -1;
            m.drugsN = 2;
            m.drugsInclude = [];
            m.sametox = false;
            m.dToxMax = 6*ones(1,size(m.dTox,2));
            m.listAllSolutions = true;
            m.changeToSubpopOnly = -1; %change pop struct to be a 100% of the nth largest subpopulation; default = 0 (keep all)
            m.popStructSet = 1:m.subpopTotal;
            m.popStructSetPrior = (1/length(m.popStructSet))*ones(size(m.popStructSet));
            m.popStructGiven = struct('N',0,'subpops',[],'contris',[]); %popstruct in abbreviated verison
            m.popFillContris = false; %if contris does not add up to 1; will use Monte Carlo to fill in the rest if popFillContris=true and montecarlo=true
            m.outputsdir = '../outputs/';
            
            %examples for m.popStructGiven
            %m.popStructGiven = m.getEmptyPopstruct(); %default empty
            %m.popStructGiven = m.readinPopstruct([outputsdir 'FILENAME']); %get popstruct from file
            %m.popStructGiven = struct('N',1,'subpops',[1],'contris',[1]); %manually define popstruct

            %derived/internal parameters
            m.drugsTotal = size(m.dEff,1);
            m.subpopTotal = size(m.dEff,2);
            
            m.results = [];
            m.resultsCombined = NaN;
            
            %not used
            m.therapyToCompare = struct('N',0,'drugs',[],'dose',[],'effects',NaN);
        end
    end
    
end


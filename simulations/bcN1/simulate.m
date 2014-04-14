% bcN1 simulation
% See README for details
% This simulation file was coded to accept input parameters to allow batch job submission
% ./runbatch.sh was run, to submit multiple jobs to cluster, running simulate_batch.m (which calls simulate.m)
% with different parameter values.
close all hidden;
clc;

%batch variables
if(~exist('subpop1pRange','var') && ~exist('changeToSubpopOnly','var'))
    error('Cannot find batch variables');
end

homedir = '../../';
matdir =  [homedir 'matlab/']; %version 2.10
datasetsdir = [matdir 'datasets/'];
matlibdir = [matdir 'lib/'];
outputsdir = './outputs/';
dataset = 'D21H30';

addpath(datasetsdir);
addpath(matlibdir);

%add cplex to path if it is not already
cplexdir = [matdir 'cplex/'];
if exist(cplexdir, 'dir'), addpath(cplexdir); end

%load data
load(dataset);

%start Modeling
m = Simulation(dsets);

%create skeleton popstructs to be filled and then optimized
%settings
subpop1pRange = str2double(subpop1pRange); %FROM BATCH VARIABLE
devRange = [1 3 5 7 9]; %difference in % between largest and second largest subpopulation
subpop1Range = 1:m.subpopTotal;
subpop2Range = 1:m.subpopTotal;

%convert % to decimal
subpop1pRange = subpop1pRange/100.0;
devRange = devRange/100.0;

%estimate popstructs size and preallocate
popstructsN = length(subpop1Range)*(length(subpop2Range)-1)*length(subpop1pRange)*length(devRange);
popstruct = m.getEmptyPopstruct();
popstructsList = repmat(popstruct,popstructsN,1);
idx = 1;

for subpop1 = subpop1Range
    for subpop2 = subpop2Range
        if(subpop1 ~= subpop2)
            for subpop1p = subpop1pRange
                for dev = devRange
                    popstruct = m.getEmptyPopstruct();
                    popstruct.N = 2;
                    popstruct.subpops = [subpop1 subpop2];
                    popstruct.contris = [subpop1p subpop1p-dev];

                    popstructsList(idx) = popstruct;
                    idx = idx + 1;
                end
            end
        end
    end
end

if(idx-1 ~= popstructsN)
    disp('Warning: Number of popstructs estimated and number enumerated do not match');
    disp(['Number of popstructs enumerated: ' num2str(idx-1) '; popstructsN: ' num2str(popstructsN) ]);
end

%settings
m.montecarlo = true;
    m.simulateN = 1;
    m.popStructSet = 1:m.subpopTotal;
    m.popStructSetPrior = []; %if empty, will default to uniform distribution
    m.subpopN = -1; %-1 means not to fix the number of subpopulations
    %m.popStructGiven = m.getEmptyPopstruct(); %default empty
    %m.popStructGiven = m.readinPopstruct([outputsdir 'FILENAME']);
    %m.popStructGiven = struct('N',1,'subpops',[1],'contris',[1]);
    m.popStructGiven = popstructsList;
    m.popFillContris = true;
m.optModel = 'efftox';
    m.drugsN = -1; %-1 means not to fix the number of drugs
    m.drugsInclude = []; %constraint to have optimal therapy to include at least the given
    m.sametox = true;
    m.dToxMax = 6*ones(1,size(m.dTox,2));
    m.listAllSolutions = false;
    m.changeToSubpopOnly = str2double(changeToSubpopOnly); %FROM BATCH VARIABLE

m.outputsdir = outputsdir;

%perform simulation
m.simulate();

%clean up
clearvars -except m;

%save workspace
save([m.outputsdir 'simulate_results' '_sametox' num2str(m.sametox) '_subpopOnly' num2str(m.changeToSubpopOnly) '_subpop1pRange' num2str(m.popStructGiven(1).contris(1)*100)], 'm');

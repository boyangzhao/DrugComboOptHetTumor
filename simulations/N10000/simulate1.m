%N10000 simulation
%See README for details
close all hidden;
clear all;
clc;

homedir = '../../';
matdir =  [homedir 'matlab/'];
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

%settings
m.montecarlo = true;
    m.simulateN = 10000;
    m.popStructSet = 1:m.subpopTotal;
    m.popStructSetPrior = []; %if empty, will default to uniform distribution
    m.subpopN = -1; %-1 means not to fix the number of subpopulations
    m.popStructGiven = m.getEmptyPopstruct(); %default empty
    %m.popStructGiven = m.readinPopstruct([outputsdir 'FILENAME']); %get popstruct from file
    %m.popStructGiven = struct('N',1,'subpops',[1],'contris',[1]); %manually define popstruct
    m.popFillContris = false;
m.optModel = 'efftox';
    m.drugsN = -1; %-1 means not to fix the number of drugs
    m.drugsInclude = []; %constraint to have optimal therapy to include at least the given
    m.sametox = true;
    m.dToxMax = 6*ones(1,size(m.dTox,2));
    m.listAllSolutions = false;
    m.changeToSubpopOnly = 0;

m.outputsdir = outputsdir;

%perform simulation
m.simulate();

%clean up
clearvars -except m;

%save workspace
save([m.outputsdir 'simulate_results' '_N' num2str(m.simulateN) '_sametox' num2str(m.sametox)], 'm');

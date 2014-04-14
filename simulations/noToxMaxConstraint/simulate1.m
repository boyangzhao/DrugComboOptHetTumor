% noToxMaxConstraint simulation
% Using given population structures
% With m.dToxMax value set high enough (e.g. 100) such that all solutions 
% along the Pareto solution will be generated
% Comment/Uncomment different m.popStructGiven for different simulations
% Change the output filename (last line of this file) accordingly
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
    m.simulateN = 1;
    m.popStructSet = 1:m.subpopTotal;
    m.popStructSetPrior = []; %if empty, will default to uniform distribution
    m.subpopN = -1; %-1 means not to fix the number of subpopulations
    
    %chk2only
    m.popStructGiven = struct('N',1,'subpops',[4],'contris',[1]);
    
    %p53only
    %m.popStructGiven = struct('N',1,'subpops',[2],'contris',[1]);
    
    %result idx 320 of N1000 (predominant subpop is p53)
    %m.popStructGiven = struct('N',14,'subpops',[2, 3, 4, 5, 6, 8, 9, 12, 13, 19, 21, 23, 25, 27],'contris',[0.55 0.02 0.01 0.01 0.01 0.15 0.09 0.01 0.05 0.01 0.04 0.01 0.03 0.01]);
    
    %result idx 2003 of N1000 (predominant subpop is chk2)
    %m.popStructGiven = struct('N',26,'subpops',[1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 25, 26, 27, 28, 30],'contris',[0.01, 0.01, 0.01, 0.55, 0.04, 0.01, 0.01, 0.01, 0.01, 0.01, 0.04, 0.01, 0.03, 0.01, 0.02, 0.01, 0.06, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.04, 0.01, 0.04]);
    
    m.popFillContris = false;
m.optModel = 'efftox';
    m.drugsN = -1; %-1 means not to fix the number of drugs
    m.drugsInclude = []; %constraint to have optimal therapy to include at least the given
    m.sametox = true;
    m.dToxMax = 100*ones(1,size(m.dTox,2));
    m.listAllSolutions = false;
    m.changeToSubpopOnly = 0;

m.outputsdir = outputsdir;

%perform simulation
m.simulate();

%clean up
clearvars -except m;

%save workspace
save([m.outputsdir 'simulate_results' '_N' num2str(m.simulateN) '_sametox' num2str(m.sametox) '_noToxMaxConstraint-chk2only'], 'm');

% Simulation with different rho values
% Using given population structure: a homogeneous shATM population with symmetric toxicity
% The rho value is not exposed for modification on purpose
% Each simulation was run with a different rho value, manually modified in OptimizationModels.m, line 100
% To also output the .lp file, uncomment line 367 in OptimizationModels.m
% Change the output filename (line 58 of this file) accordingly
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
m.montecarlo = false;
    m.simulateN = 1;
    m.popStructSet = 1:m.subpopTotal;
    m.popStructSetPrior = []; %if empty, will default to uniform distribution
    m.subpopN = -1; %-1 means not to fix the number of subpopulations
    %m.popStructGiven = m.getEmptyPopstruct(); %default empty
    %m.popStructGiven = m.readinPopstruct([outputsdir 'FILENAME']); %get popstruct from file
    m.popStructGiven = struct('N',1,'subpops',[3],'contris',[1]); %manually define popstruct
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
save([m.outputsdir '10e-0'], 'm');

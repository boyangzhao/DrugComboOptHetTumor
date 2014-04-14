%Analyzes noToxMaxConstraint simulation results 
close all hidden;
clear all;
clc;

homedir = '../../';
matlibdir = [homedir 'matlab/lib/'];

addpath(matlibdir);

%% Analysis parameters
%Uncomment/comment files below to analyze select simulation results
load './simulate_results_N1_sametox1_noToxMaxConstraint-result2003';
%load './simulate_results_N1_sametox1_noToxMaxConstraint-result320';
%load './simulate_results_N1_sametox1_noToxMaxConstraint-p53only';
%load './simulate_results_N1_sametox1_noToxMaxConstraint-chk2only';

%% Analyze simulation results
m.generatePlotsPerSol(m.results(1));

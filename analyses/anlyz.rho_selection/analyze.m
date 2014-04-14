%Analyze rho_selection simulation results
close all hidden;
clear all;
clc;

homedir = '../../';
matlibdir = [homedir 'matlab/lib/'];

addpath(matlibdir);

%% Analysis parameters
%change filename (e.g. ./10e-9) to analyze simulation result of interest
load './10e-9';

%% Analyze simulation results
m.generatePlotsPerSol(m.results(1));

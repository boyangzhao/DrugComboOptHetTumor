%Analyzes N10000_sametox0 simulation results
%Per solution analysis
close all hidden;
clear all;
clc;

homedir = '../../';
matlibdir = [homedir 'matlab/lib/'];
outputsdir = '../outputs/';

addpath(matlibdir);
addpath(outputsdir);

load 'simulate_results_N10000_sametox0';

%Before run: change sortsols = true (in Analysis.m)
m.generatePlotsPerSol(m.results(320));
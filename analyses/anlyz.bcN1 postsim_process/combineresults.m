%Combines simulation results from bcN1
close all hidden;
clear all;
clc;

homedir = '../../';
matlibdir = [homedir 'matlab/lib/'];
outputsdir = '../outputs_bc/';

addpath(matlibdir);
addpath(outputsdir);

load 'simulate_results_sametox1_subpopOnly0_subpop1pRange50'; m1=m;
load 'simulate_results_sametox1_subpopOnly0_subpop1pRange45'; m2=m;
load 'simulate_results_sametox1_subpopOnly0_subpop1pRange40'; m3=m;
m = Analysis.combineResults([m1;m2;m3], 'subpop1pRange50+subpop1pRange40');
save(['./' 'simulate_results_sametox1_subpopOnly0_combined'], 'm');
clearvars m m1 m2 m3 m_combined

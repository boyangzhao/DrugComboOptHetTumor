%Compares between N10000 and (N10000_subpopOnly1 & N10000_subpopOnly2) simulation results
close all hidden;
clear all;
clc;

homedir = '../../';
matlibdir = [homedir 'matlab/lib/'];
outputsdir = '../outputs/';

addpath(matlibdir);
addpath(outputsdir);

load 'simulate_results_N10000_sametox1_processed'; m0=m;
load 'simulate_results_N10000_sametox1_subpopOnly1'; m1=m;
load 'simulate_results_N10000_sametox1_subpopOnly2'; m2=m;

[consensus,~,~]=Analysis.getConsensus(m0,[m1;m2],true);

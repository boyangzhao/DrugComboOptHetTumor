%Analysis of N10000 simulation results
%Global statistical and sensitivity analyses
close all hidden;
clear all;
clc;

homedir = '../../';
matlibdir = [homedir 'matlab/lib/'];
outputsdir = '../outputs/';

addpath(matlibdir);
addpath(outputsdir);

load 'simulate_results_N10000_sametox1_processed';

%Change 'highest' to other values accordingly - to examine other regimens
m.analyze('highest'); %lowest, highest, compromise, or index value
m.generatePlots();
%Post simulation processing - removes any weak pareto optimal solutions, if exists
close all hidden;
clear all;
clc;

homedir = '../../';
matfilesdir = [homedir 'matlab/'];
outputsdir = '../outputs/';

addpath(matfilesdir);

load([outputsdir 'simulate_results_N10000_sametox1']);
weakParetoSet = m.removeWeakParetoOptimal();

disp('**********************')
disp(['Number of weak pareto solutions: ' num2str(length(weakParetoSet))]);
weakParetoSet
disp(' ');
if(~isempty(weakParetoSet))
    save(['./' 'simulate_results_N10000_sametox1_processed'], 'm');
end

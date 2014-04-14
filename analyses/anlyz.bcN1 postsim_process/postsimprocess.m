%Post simulation processing - removes any weak pareto optimal solutions, if exists
close all hidden;
clear all;
clc;

homedir = '../../';
matfilesdir = [homedir 'matlab/'];
outputsdir = '../outputs_bc/';

addpath(matfilesdir);

for subpopOnly = 0
    load([outputsdir 'simulate_results_sametox1_subpopOnly' num2str(subpopOnly)  '_combined']);
    weakParetoSet = m.removeWeakParetoOptimal();
    
    disp('**********************')
    disp(['Number of weak pareto solutions for subpopOnly' num2str(subpopOnly) ': ' num2str(length(weakParetoSet))]);
    weakParetoSet
    disp(' ');
    if(~isempty(weakParetoSet))
        save(['./' 'simulate_results_sametox1_subpopOnly' num2str(subpopOnly) '_combined_processed'], 'm');
    end
end

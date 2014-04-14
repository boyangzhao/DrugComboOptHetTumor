%Pairwise correlation comparison of drug frequencies among different regimens on Pareto frontier
close all hidden;
clear all;
clc;

homedir = '../../';
matlibdir = [homedir 'matlab/lib/'];
outputsdir = '../outputs/';

addpath(matlibdir);
addpath(outputsdir);

freqs = zeros(21,6);
for i = 2:7
    load 'simulate_results_N10000_sametox1_processed';
    m.analyze(i);
    freqs(:,i-1) = sum(m.analysis.corr_drug,1)';    
    clear m;
end

freq_corr = zeros(6,6);
for i = 1:6
    for j = 1:6
        freq_corr(i,j) = corr(freqs(:,i),freqs(:,j),'type','spearman');
    end
end

t = {'1';'2';'3';'4';'5';'6'}
Utilities.heatmap_m(freq_corr, ...
                    'rowlabels', t,...
                    'collabels', t,...
                    'xlabel', 'Solution',...
                    'ylabel', 'Solution',...
                    'range',[0.5 1],...
                    'cmap',Utilities.getcmap('hot'),...
                    'cbar',true);
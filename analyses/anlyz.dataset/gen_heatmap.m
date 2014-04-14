close all hidden;
clear all;
clc;

homedir = '../../';
matdsdir = [homedir 'matlab/datasets/'];
addpath(matdsdir);
load D21H30;

%Generates clustergram of D21H30 efficacy dataset
clustergram(dEff, ...
          'rowlabels',drug_labels_long, ...
          'columnlabels',subpop_labels, ...
          'standardize',3, ...
          'colormap',getcmap('rwb'), ...
          'displayrange', 6);

%Generates heatpmap of D21H30 toxicity dataset
heatmap_m(dTox, ...
          'rowlabels',drug_labels_long, ...
          'collabels',dTox_labels, ...
          'cbarlabel','Toxicity', ...
          'rotateXLabel', 45, ...
          'cmap',getcmap('gray_binary'), ...
          'cbar', false);
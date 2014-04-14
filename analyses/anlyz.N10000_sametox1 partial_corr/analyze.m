%Matthews correlation and partial correlation analysis of N10000, N10000_subpopOnly1 results
close all hidden;
clear all;
clc;

homedir = '../../';
matlibdir = [homedir 'matlab/lib/'];
outputsdir = '../outputs/';
addpath(matlibdir);
addpath(outputsdir);

%% Analysis parameters
%Change parameter value here for choosing which regiment to analyze
drugsNcombo = 6;

%Comment/uncomment the following simulation results for analyzing either N10000, N10000_subpopOnly1
%load 'simulate_results_N10000_sametox1_processed';
load 'simulate_results_N10000_sametox1_subpopOnly1';

%% Derived parameters
regimen = drugsNcombo+1; %regimen also includes a 0-drug combo; so regimen and drugsNcombo is offset by 1
m.analyze(regimen);

%% Partial correlation analysis
dmatrix = m.analysis.popNdrug_matrix(m.analysis.drugsToInclude,:)'; %row = subpopN; column = drug
labels = m.analysis.drugsIncluded_labels;

[rho,pval] = corr(dmatrix, 'type', 'spearman');
[part_rho,part_pval] = partialcorr(dmatrix, [1:m.subpopTotal]', 'type', 'spearman');

% calculate critical rho value
N = size(dmatrix,1);
df = N-2;
rcritical = Utilities.calcRcritical(0.05, df);
disp('Drug freq vs tumor complexity');
disp(['DF: ' num2str(df)]);
disp(['Critical rho value at p < 0.05: ' num2str(rcritical)]);
disp(' ');

c = clustergram(rho,...
            'RowLabels',labels, ...
            'ColumnLabels',labels, ...
            'colormap',Utilities.getcmap('rwb'),...
            'standardize',3, ...
            'ShowDendrogram', 'on',...
            'linkage','average',...
            'rowpdist','euclidean',...
            'columnpdist','euclidean');
addTitle(c, 'Correlation matrix of drugs based on tumor complexity', 'FontSize', 12);

c = clustergram(part_rho,...
            'RowLabels',labels, ...
            'ColumnLabels',labels, ...
            'colormap',Utilities.getcmap('rwb'),...
            'standardize',3,...
            'ShowDendrogram', 'on',...
            'linkage','average',...
            'rowpdist','euclidean',...
            'columnpdist','euclidean');
addTitle(c, 'Partial correlation matrix of drugs, controlled for tumor complexity', 'FontSize', 12);

%get upper triangle matrix, excluding diagonal; use this as logical matrix
%to select the unique correlation values
%{
select = logical(triu(ones(size(rho)),1));
rho_vals = rho(select);
part_rho_vals = part_rho(select);
%}

%% Matthews correlation
resultsN = size(m.results,1);
drugsMatrix = zeros(resultsN, m.drugsTotal); %binary matrix; row = population; column = different drugs; 
for i = 1:resultsN
    drugsMatrix(i,:) = m.results(i).sols(regimen,:);
end

drugsMatrix = drugsMatrix(:,m.analysis.drugsToInclude);
[drho,dpval] = corr(drugsMatrix, 'type', 'pearson'); %pearson for dichotomous data, reduced to Matthews correlation

c = clustergram(drho,...
            'RowLabels',labels, ...
            'ColumnLabels',labels, ...
            'colormap',Utilities.getcmap('rwb'),...
            'standardize',3, ...
            'ShowDendrogram', 'on',...
            'linkage','average',...
            'rowpdist','euclidean',...
            'columnpdist','euclidean');
addTitle(c, 'Matthews correlation', 'FontSize', 12);

%% Calculate critical rho values
N = size(drugsMatrix,1);
df = N-2;
rcritical = Utilities.calcRcritical(0.05, df);
disp('Drug choice vs tumor populations');
disp(['DF: ' num2str(df)]);
disp(['Critical rho value at p < 0.05: ' num2str(rcritical)]);

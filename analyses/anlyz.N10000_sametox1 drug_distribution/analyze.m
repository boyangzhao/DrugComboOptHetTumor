%Analysis of N10000,  N10000_subpopOnly1, and N10000_subpopOnly2 results
%decomposed based on predominant and second largest subpopulation proportions
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
drugsNcombo = 6; %number of drugs in drug combination

%Comment/uncomment files to analyze different simulation results;
%Change inputs to Analysis.getConsensus accordingly
load 'simulate_results_N10000_sametox1_processed'; m0=m;
load 'simulate_results_N10000_sametox1_subpopOnly1'; m1=m;
load 'simulate_results_N10000_sametox1_subpopOnly2'; m2=m;

[consensus,~,~]=Analysis.getConsensus(m0,[m1;m2]);

%% Get distribution info of subpopulations
resultsN = size(m0.results,1); %max, second max, third max
sub_percentages = zeros(resultsN,3);
fmatrix = zeros(101,101); %frequency matrix; row = predominant subpopulation %s (0-100% (index 1 to 101)); col = deviation % (second largest subpop - largest subpop %s)
cmatrix = zeros(101,101); %consensus matrix; frequency of pop that have different drug combinations
cpmatrix = zeros(101,101); %consensus matrix compared to single best drug of predominant; fraction of pop not containing single-best drug for predominant

for i = 1:resultsN
    [v,idx] = sort(m0.results(i).pop);
    sub_percentages(i,:) = [v(end), v(end-1), v(end-2)];
    
    subpop1 = uint8(v(end)*100);
    subpop2= uint8(v(end-1)*100);
    
    row = subpop2+1;
    col = subpop1+1;
    fmatrix(row,col) = fmatrix(row,col) + 1;%add to frequency deviation % of current population
    
    if(~consensus(i,drugsNcombo)) %drug combo differ b/t optimization based on entire heterogeneity versus just predominant
        sol = m1.results(i).sols(2,:); %solution for 1-drug combination
        drug_idx = find(sol == 1);
        if(m0.results(i).sols(drugsNcombo+1,drug_idx) == 0)
            %sol does not contain single-best drug for predominant
            cpmatrix(row,col) = cpmatrix(row,col) + 1;
        end
        cmatrix(row,col) = cmatrix(row,col) + 1;
    end
end
cmatrix_normalized = cmatrix./fmatrix; %normalize to frequency of pops in each block
cpmatrix_normalized = cpmatrix./cmatrix; %calculate fraction of pop (with different drug combo) that does not contain single-best for predominant
cpmatrix_normalized(isnan(cpmatrix_normalized)) = 0;

%% Plot freq distributions of subpopulation %s
figure;
hist(sub_percentages(:,1),0.025:0.05:1);
set(gca,'xlim',[0 1]);
title('Predominant subpopulation');
xlabel('Predominant subpopulation %');
ylabel('Frequency');

figure;
hist(sub_percentages(:,2),0.025:0.05:1);
set(gca,'xlim',[0 1]);
title('Second highest subpopulation');
xlabel('Second largest subpopulation %');
ylabel('Frequency');

%flip matrice when plotting heatmap so for y axis, it goes from 0 (bottom) to 100 (top) instead of 100 to 0
%heatmap pop freq; plotted as predominant (largest) subpopulation %s vs. second largest subpopulation %s
figure;
Utilities.heatmap_m(flipud(fmatrix),...
                              'ylabel','Second largest subpopulation %',...
                              'xlabel','Predominant subpopulation %',...
                              'collabels',0:100,...
                              'rowlabels',100:-1:0,...
                              'range',[0 30],...
                              'cmap',Utilities.getcmap('hot',true));

%heatmap frac not containing single-best drug for predominant subpop; plotted as predominant (largest) subpopulation %s vs. second largest subpopulation %s
%get frac of single-best drug for predominant subpop
figure;
Utilities.heatmap_m(flipud(cmatrix_normalized),...
                              'ylabel','Second largest subpopulation %',...
                              'xlabel','Predominant subpopulation %',...
                              'collabels',0:100,...
                              'rowlabels',100:-1:0,...
                              'range',[0 1],...
                              'cmap',Utilities.getcmap('rwb',false,[0.5 1]),...
                              'title','cmatrix-normalized');

figure;
Utilities.heatmap_m(flipud(cpmatrix_normalized),...
                              'ylabel','Second largest subpopulation %',...
                              'xlabel','Predominant subpopulation %',...
                              'collabels',0:100,...
                              'rowlabels',100:-1:0,...
                              'range',[0 1],...
                              'cmap',Utilities.getcmap('rwb',false,[0.5 1]),...
                              'title','cpmatrix-normalized');

%% Analyze results along the two axes: 1) predominant %, versus 2) deviation % from predominant %)
%initialize matrices now rotated with row = predominant % (0-100%); col = deviation % from predominant % (0-100%)
%initialize first with -1s, later change all the -1s to be NaNs
cmatrix_norm_rot = -1*ones(101,101);
cpmatrix_norm_rot = -1*ones(101,101);
for i = 0:100 %referring to %
    for j = 0:100
        %coord tracks index in the matrix, which is +1 of the actual %;
        %since starting index in matrix/vector is 1 and not 0, but percentages start at 0
        coord =struct('subpop1',i+1, 'subpop2', j+1,'dev',i-j+1); %subpop1 = predominant %; subpop2 = second largest subpop %
        dev = i-j;
        if(dev >= 0 && fmatrix(coord.subpop2, coord.subpop1) > 1)
            cmatrix_norm_rot(coord.subpop1, coord.dev) = cmatrix_normalized(coord.subpop2, coord.subpop1);
            cpmatrix_norm_rot(coord.subpop1, coord.dev) = cpmatrix_normalized(coord.subpop2, coord.subpop1);
        end
    end
end

cmatrix_norm_rot(cmatrix_norm_rot<0) = NaN;
cpmatrix_norm_rot(cpmatrix_norm_rot<0) = NaN;

%calculate variance for each group first; can still result in NaNs if all elements in group are NaNs
SSc_subpop1 = nanvar(cmatrix_norm_rot,0,1);
SSc_dev = nanvar(cmatrix_norm_rot,0,2);
SScp_subpop1 = nanvar(cpmatrix_norm_rot,0,1);
SScp_dev = nanvar(cpmatrix_norm_rot,0,2);

figure;
Utilities.heatmap_m(flipud(cmatrix_norm_rot),...
                              'ylabel','Predominant subpopulation %',...
                              'xlabel','Difference (in %) between predominant and second largest subpopulations',...
                              'collabels',0:100,...
                              'rowlabels',100:-1:0,...
                              'range',[0 1],...
                              'cmap',Utilities.getcmap('rwb',false,[0.5 1]),...
                              'title','cmatrix-norm-rot');
                          
figure;
Utilities.heatmap_m(flipud(cpmatrix_norm_rot),...
                              'ylabel','Predominant subpopulation %',...
                              'xlabel','Difference (in %) between predominant and second largest subpopulations',...
                              'collabels',0:100,...
                              'rowlabels',100:-1:0,...
                              'range',[0 1],...
                              'cmap',Utilities.getcmap('rwb',false,[0.5 1]),...
                              'title','cpmatrix-norm-rot');

%Analysis of bcN1,  bcN1_subpopOnly1, and bcN1_subpopOnly2 results
%decomposed based on predominant and second largest subpopulation proportions
close all hidden;
clear all;
clc;

homedir = '../../';
matlibdir = [homedir 'matlab/lib/'];
outputsdir = '../outputs_bc/';

addpath(matlibdir);
addpath(outputsdir);

%% Analysis parameters (can change)
%Comment/uncomment files to analyze different simulation results;
%Change inputs to Analysis.getConsensus accordingly
load 'simulate_results_sametox1_subpopOnly0_combined_processed'; m0=m;
load 'simulate_results_sametox1_subpopOnly1'; m1=m;
load 'simulate_results_sametox1_subpopOnly2'; m2=m;

[consensus,~,~]=Analysis.getConsensus(m0,[m1;m2]);

%% Analysis parameters
drugsNcombos = 1:6; %number of drugs in drug combination

%% Get distribution info of subpopulations
resultsN = size(m0.results,1); %max, second max, third max
sub_percentages = zeros(resultsN,3);
fmatrix = zeros(101,101); %frequency matrix; row = predominant subpopulation %s (0-100% (index 1 to 101)); col = deviation % (second largest subpop - largest subpop %s)

cmatrix{1} = zeros(101,101); %consensus matrix; frequency of pop that have different drug combinations
cpmatrix{1} = zeros(101,101); %consensus matrix compared to single best drug of predominant; fraction of pop not containing single-best drug for predominant
cmatrix = repmat(cmatrix,length(drugsNcombos),1); %cell array, each cell is a different regimen
cpmatrix = repmat(cpmatrix,length(drugsNcombos),1); %cell array, each cell is a different regimen

for i = 1:resultsN
    [v,idx] = sort(m0.results(i).pop);
    sub_percentages(i,:) = [v(end), v(end-1), v(end-2)];
    
    subpop1 = uint8(v(end)*100);
    subpop2= uint8(v(end-1)*100);
    
    row = subpop2+1;
    col = subpop1+1;
    fmatrix(row,col) = fmatrix(row,col) + 1;%add to frequency deviation % of current population
    
    for n = 1:length(drugsNcombos)
        if(~consensus(i,drugsNcombos(n))) %drug combo differ b/t optimization based on entire heterogeneity versus just predominant
            sol = m1.results(i).sols(2,:); %solution for 1-drug combination
            drug_idx = find(sol == 1);
            if(m0.results(i).sols(drugsNcombos(n)+1,drug_idx) == 0)
                %sol does not contain single-best drug for predominant
                cpmatrix{drugsNcombos(n)}(row,col) = cpmatrix{drugsNcombos(n)}(row,col) + 1;
            end
            cmatrix{drugsNcombos(n)}(row,col) = cmatrix{drugsNcombos(n)}(row,col) + 1;
        end
    end
end

fmatrix_rep{1} = fmatrix;
fmatrix_rep = repmat(fmatrix_rep,length(drugsNcombos),1);
cmatrix_normalized = cellfun(@(x,y) x./y, cmatrix, fmatrix_rep, 'UniformOutput', false);
cpmatrix_normalized = cellfun(@(x,y) x./y, cpmatrix, cmatrix, 'UniformOutput', false);
%{
for n = 1:length(drugsNcombos) %set the NaNs to zero; NaNs are coming from divisions by 0 in operations above
    cmatrix_normalized{drugsNcombos(n)}(isnan(cmatrix_normalized{drugsNcombos(n)})) = 0;
    cpmatrix_normalized{drugsNcombos(n)}(isnan(cpmatrix_normalized{drugsNcombos(n)})) = 0;
end
%}

%% Plot freq distributions of subpopulation %s
%flip matrice when plotting heatmap so for y axis, it goes from 0 (bottom) to 100 (top) instead of 100 to 0
%heatmap pop freq; plotted as predominant (largest) subpopulation %s vs. second largest subpopulation %s
figure(1);
Utilities.heatmap_m(flipud(fmatrix),...
                              'ylabel','Second largest subpopulation %',...
                              'xlabel','Predominant subpopulation %',...
                              'collabels',0:100,...
                              'rowlabels',100:-1:0,...
                              'range',[0 30],...
                              'cmap',Utilities.getcmap('hot',true));


%% Generate/plot matrices of predominant % versus difference % (between predominant and second largest subpopulation)
%initialize matrices now rotated with row = predominant % (0-100%); col = deviation % from predominant % (0-100%)
%initialize first with -1s, later change all the -1s to be NaNs
fmatrix_rot = zeros(101,101);
cmatrix_normalized_rot{1} = fmatrix_rot;
cpmatrix_normalized_rot{1} = fmatrix_rot;
cmatrix_normalized_rot = repmat(cmatrix_normalized_rot,length(drugsNcombos),1);
cpmatrix_normalized_rot = repmat(cpmatrix_normalized_rot,length(drugsNcombos),1);

for i = 0:100 %referring to %
    for j = 0:100
        %coord tracks index in the matrix, which is +1 of the actual %;
        %since starting index in matrix/vector is 1 and not 0, but percentages start at 0
        coord =struct('subpop1',i+1, 'subpop2', j+1,'dev',i-j+1); %subpop1 = predominant %; subpop2 = second largest subpop %
        dev = i-j;
        if(dev >= 0 && fmatrix(coord.subpop2, coord.subpop1) > 1)
            fmatrix_rot(coord.subpop1, coord.dev) = fmatrix(coord.subpop2, coord.subpop1);
            
            for n = 1:length(drugsNcombos)
                cmatrix_normalized_rot{drugsNcombos(n)}(coord.subpop1, coord.dev) = cmatrix_normalized{drugsNcombos(n)}(coord.subpop2, coord.subpop1);
                cpmatrix_normalized_rot{drugsNcombos(n)}(coord.subpop1, coord.dev) = cpmatrix_normalized{drugsNcombos(n)}(coord.subpop2, coord.subpop1);
            end
        end
    end
end

figure(2);
Utilities.heatmap_m(flipud(fmatrix_rot),...
                              'ylabel','Predominant subpopulation %',...
                              'xlabel','Difference (in %) between predominant and second largest subpopulations',...
                              'collabels',0:100,...
                              'rowlabels',100:-1:0,...
                              'range',[0 1000],...
                              'cmap',Utilities.getcmap('rwb',false,[0.5 1]),...
                              'title','fmatrix-rot');
axis([0.5 11.5 49.5 62.5]);

%% Go through each element in the fmatrix_rot, and generate consensus if frequency > 1
subplotnum = 1;
for subpop1 = 1:size(fmatrix_rot,1) %subpop1 index
    for dev = 1:size(fmatrix_rot,2) %dev index; actual dev % is one less
        if(fmatrix_rot(subpop1,dev) > 1)
            c_fraction = zeros(length(drugsNcombos),2);
            cp_fraction = zeros(length(drugsNcombos),2);
            for n = 1:length(drugsNcombos)
                c_fraction(n,1) = cmatrix_normalized_rot{drugsNcombos(n)}(subpop1,dev);
                cp_fraction(n,1) = cpmatrix_normalized_rot{drugsNcombos(n)}(subpop1,dev);
            end
            
            c_fraction(:,2) = 1 - c_fraction(:,1);
            cp_fraction(:,2) = 1 - cp_fraction(:,1);
            
            figure(3);
            subplot(3,5,subplotnum);
            h = bar(1:length(drugsNcombos), c_fraction, 1, 'stack');
            set(gca,'XLim',[0.5 length(drugsNcombos)+0.5], ...
                        'XTick',[1:length(drugsNcombos)], ...
                        'XTickLabel',1:length(drugsNcombos));
            %legend('Different','Same');
            title(['c-fraction - ' 'subpop1: ' num2str(subpop1) '; dev: ' num2str(dev-1)]);
            xlabel('Regimens');
            ylabel('Frequency');
            Utilities.applyPlotBasicStyles(gcf);
            set(h(1),'Facecolor',[59/255 158/255 59/255]); %different combo
            set(h(2),'Facecolor',[50/255 108/255 191/255]); %same combo
            
            figure(4);
            subplot(3,5,subplotnum);
            h = bar(1:length(drugsNcombos), cp_fraction, 1, 'stack');
            set(gca,'XLim',[0.5 length(drugsNcombos)+0.5], ...
                        'XTick',[1:length(drugsNcombos)], ...
                        'XTickLabel',1:length(drugsNcombos));
            %legend('Different','Same');
            title(['cp-fraction - ' 'subpop1: ' num2str(subpop1) '; dev: ' num2str(dev-1)]);
            xlabel('Regimens');
            ylabel('Frequency');
            Utilities.applyPlotBasicStyles(gcf);
            set(h(1),'Facecolor',[50/255 50/255 50/255]); %does not contain single-best
            set(h(2),'Facecolor',[200/255 200/255 200/255]); %contains single-best
            
            subplotnum = subplotnum+1;
    end
    end
end

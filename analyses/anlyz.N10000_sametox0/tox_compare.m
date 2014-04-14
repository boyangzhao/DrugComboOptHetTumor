%Compares N10000_sametox0 and N10000 simulation results
close all hidden;
clear all;
clc;

homedir = '../../';
matlibdir = [homedir 'matlab/lib/'];
outputsdir = '../outputs/';

addpath(matlibdir);
addpath(outputsdir);

%% compare sametox0 and sametox1
% tally the toxicity profile for each regimen of each pop struct
load 'simulate_results_N10000_sametox0';
dTox = m.dTox;
toxprofile_sametox0 = getToxProfile(dTox, m.results);
clearvars m;

load 'simulate_results_N10000_sametox1_processed';
toxprofile_sametox1 = getToxProfile(dTox, m.results);
clearvars m;

%% normal distr fit on sametox0 and sametox1
load 'simulate_results_N10000_sametox0';
for i = 1:6 %loop for drugsN
    for j = 1:length(m.dTox_labels)
        [mu sigma] = normfit(toxprofile_sametox1(i).tp(:,j));
        toxprofile_sametox1(i).normfit = [toxprofile_sametox1(i).normfit;mu sigma];
    end
end

for i = 1:6 %loop for drugsN
    for j = 1:length(m.dTox_labels)
        [mu sigma] = normfit(toxprofile_sametox0(i).tp(:,j));
        toxprofile_sametox0(i).normfit = [toxprofile_sametox0(i).normfit;mu sigma];
    end
end

%% plot sametox0 and sametox1
toxIndices = [1,6]; %toxicties to analyze; each value = index in toxicity profile matrix, same order as dTox (column index)

for toxIdx = toxIndices
    figure;
    for i = 1:6
        subplot(2,3,i)
        [n x] = hist(toxprofile_sametox1(i).tp(:,toxIdx),0:6);
        bar(x, n/sum(n),1);
        set(gca,'XLim',[-0.5 6.5]);
        set(gca,'YLim',[0 1]);
        title(['same toxicity; drugN = ' num2str(i)]);
    end

    figure;
    for i = 1:6
        subplot(2,3,i)
        [n x] = hist(toxprofile_sametox0(i).tp(:,toxIdx),0:6);
        bar(x, n/sum(n),1);
        set(gca,'XLim',[-0.5 6.5]);
        set(gca,'YLim',[0 1]);
        title(['actual toxicity; drugN = ' num2str(i)]);
    end

    figure;
    hold on;
    for i = 1:6
        x = 0:0.1:6;
        h = plot3(i*ones(length(x),1), x, normpdf(x, toxprofile_sametox0(i).normfit(toxIdx,1), toxprofile_sametox0(i).normfit(toxIdx,2)), 'r', ...
                  i*ones(length(x),1), x, normpdf(x, toxprofile_sametox1(i).normfit(toxIdx,1), toxprofile_sametox1(i).normfit(toxIdx,2)), 'k');
        set(h,'LineWidth',2);
    end
    hold off;
    xlabel('No. of drugs');
    ylabel('No. of overlaps in toxicity');
    zlabel('Relative frequency');
    set(gca,'XTick',1:6);
    set(gca,'YTick',0:6);
    set(gca,'ZTick',0:1);
    set(gca,'XLim',[1 6]);
    set(gca,'YLim',[0 6]);
    set(gca,'ZLim',[0 1]);
    set(gca,'XGrid','on','YGrid','on','ZGrid','off');
    Utilities.applyPlotBasicStyles(gcf);
    title(m.dTox_labels(toxIdx));
    legend(h,{'Different toxicity','Same toxicity'});
    view([-180 -100 450]);
end

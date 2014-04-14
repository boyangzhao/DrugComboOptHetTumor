%% Optimization Models Class for combination therapy optimization on heterogeneous populations
% Author: Boyang Zhao
% Broad analyses of simulation results
% Methods associated with an instance perform broad analyses of the given instance
% Static methods analyzes across instances; require multiple instances as inputs
% Dependencies: corr_pb, heatmap_m

classdef Analysis < handle & OptimizationModels

    properties (SetAccess = private, GetAccess = private)
        figN;
    end
    
    properties (SetAccess = private, GetAccess = public)
        analysis; %analyzed results
    end
    
    methods (Access = public)
        function m = Analysis(dsets)
            %Constructor
            m = m@OptimizationModels(dsets);
            m.figN = 0;
            m.analysis = [];
        end
        
        function analyze(m, solIdx)
            %analysis of result (current instance of m)
            resultsN = size(m.results,1);
            
            m.analysis = struct('drugsN', zeros(1,m.drugsTotal+1), ...
                                'drugs', zeros(1,m.drugsTotal), ...
                                'subpopN', zeros(1,m.subpopTotal), ...
                                'subpop', zeros(1,m.subpopTotal), ...
                                ...
                                'corr_pop', -2*ones(resultsN, m.subpopTotal), ...
                                'corr_drug', -2*ones(resultsN, m.drugsTotal), ...
                                'popNdrug_matrix', zeros(m.drugsTotal, m.subpopTotal), ...
                                'corr_popNdrugN', zeros(resultsN, 2), ...
                                'corr_popNdrug', zeros(1,m.drugsTotal), ...
                                ...
                                'drugsToInclude', logical(ones(1,m.drugsTotal)), ... %default to include all drugs in analysis
                                'drugsIncluded_labels', [], ...
                                'corr', -2*ones(m.drugsTotal, m.subpopTotal), ...
                                'corr_variables', -2*ones(m.drugsTotal,5), ...
                                'corr_variables_labels', [], ...
                                'corr_matrix', -2*ones(5,5));
                            
            %To avoid issues with NaNs in correlations, any drugs with
            %frequencies equal to 0 are removed (using an inclusion
            %variable m.analysis.drugsToInclude). This affects the
            %following variables: corr, corr_variables, and corr_matrix
            
            for i = 1:resultsN
                %assumes: m.getSol(i,solIdx) returns a vector
                
                %distribution of number of drugs
                c = sum(m.getSol(i,solIdx),2);
                n = hist(c,0:m.drugsTotal);
                m.analysis.drugsN = m.analysis.drugsN + n;
                
                %distribution of drugs used
                c = sum(m.getSol(i,solIdx),1);
                m.analysis.drugs = m.analysis.drugs + c;
                
                %distribution of number of subpopulations
                c = sum(m.results(i).pop>0,1);
                m.analysis.subpopN(c) = m.analysis.subpopN(c) + 1;
                
                %popN drugN matrix (each row = drug; each column = number of subpopulations)
                m.analysis.popNdrug_matrix(:,c) = m.analysis.popNdrug_matrix(:,c) + m.getSol(i,solIdx)';
                
                %distribution of subpopulation used
                m.analysis.subpop = m.analysis.subpop + (m.results(i).pop>0)';
                
                %for correlation; row = result; column = drug/pop
                m.analysis.corr_pop(i,:) = m.results(i).pop';
                m.analysis.corr_drug(i,:) = m.getSol(i,solIdx);
            end
            
            %correlation between popN and drugN
            i = sum(m.analysis.corr_pop>0,2);
            j = sum(m.analysis.corr_drug,2);
            m.analysis.corr_popNdrugN = [i j];
            
            %correlation between popN and drug
            m.analysis.corr_popNdrug = corr((1:m.subpopTotal)', ...
                                            m.analysis.popNdrug_matrix', ...
                                            'type','spearman');
                                        
            %get indices of drugs with freq > 0; drugs with freq = 0 will be excluded
            m.analysis.drugsToInclude = (sum(m.analysis.corr_drug,1) > 0);
            
            m.analysis.drugsIncluded_labels = m.drug_labels(m.analysis.drugsToInclude);
            
            %sensitivity analysis
            for i = 1:m.drugsTotal
                if(m.analysis.drugsToInclude(i))
                    for j = 1:m.subpopTotal
                        m.analysis.corr(i,j) = Utilities.corr_pb(m.analysis.corr_drug(:,i), m.analysis.corr_pop(:,j));
                    end
                end
            end
            
            m.analysis.corr =  m.analysis.corr(m.analysis.drugsToInclude,:);
            m.analysis.corr_variables = -2*ones(sum(m.analysis.drugsToInclude),5);
            
            m.analysis.corr_variables_labels = {'avg corr';'kurtosis corr';'freq';'avg dEff';'range corr'};
            m.analysis.corr_variables(:,1) = mean(m.analysis.corr,2);
            m.analysis.corr_variables(:,2) = kurtosis(m.analysis.corr,0,2);
            m.analysis.corr_variables(:,3) = sum(m.analysis.corr_drug(:,m.analysis.drugsToInclude),1)';
            m.analysis.corr_variables(:,4) = mean(m.dEff(m.analysis.drugsToInclude,:),2);
            m.analysis.corr_variables(:,5) = range(m.analysis.corr,2);
            
            m.analysis.corr_matrix = corr(m.analysis.corr_variables,'type','spearman');
        end
        
        function weakParetoSet = checkWeakParetoOptimal(m)
            %check for weak Pareto optimal solutions in results
            %output list of indices (in m.results) that contains weak Pareto optimal solutions
            weakParetoSet = [];
            for idx = 1:length(m.results)
                objN = size(m.results(idx).objvals,2);
                for obj = 1:objN
                    num_total = length(m.results(idx).objvals(:,obj));
                    num_unique = length(unique(m.results(idx).objvals(:,obj)));
                    
                    if(~Utilities.comprval(num_total, num_unique))
                        if(num_unique < num_total)
                            weakParetoSet = [weakParetoSet;idx];
                        else
                            disp('Warning: Something is wrong, total number of solutions is less number of unique solutions');
                        end
                    end
                end
            end
        end
        
        function weakParetoSet = removeWeakParetoOptimal(m)
            %REQUIRES: 1) the objective vals are sorted - this will be the case based the optimization approach used in OptimizationModels
            %                       objective vals are in ascending order with increasing number of drugs in the combo.
            %                       This is required so for each redundant value (containing weak Pareto optimal solutions), the last solution 
            %                       would be Pareto optimal (if the obj vals are sorted).
            %DEPENDENCIES: calls Analysis.checkWeakParetoOptimal to get the list of weak Pareto optimal solutions
            %MODIFIES: m.results.objvals and m.results.sols
            weakParetoSet = m.checkWeakParetoOptimal();
            
            for idx = 1:length(weakParetoSet)
                resultIdx = weakParetoSet(idx);
                objN = size(m.results(resultIdx).objvals,2);
                for obj = 1:objN
                    idx_all = 1:length(m.results(resultIdx).objvals(:,obj));
                    [~,idx_unique,~] = unique(m.results(resultIdx).objvals(:,obj));
                    
                    %the 'ia' output from unique function doesn't seem to be stable, the index chosen depending on whether setOrder is 
                    %explicitly defined as 'sorted' or not will result in different indices. Below method will manually find the indices of 
                    %redundant values and only keep the last index per redunant value
                    idx_redundant = setxor(idx_all,idx_unique); %list of indices (in objvals) containing the redundant value
                    if(~isempty(idx_redundant)) %otherwise there are no redundant values for this objective
                        redundantVals = m.results(resultIdx).objvals(idx_redundant, obj); %list of redundant values (in objvals)
                        IdxToExclude = [];
                        for r = 1:length(redundantVals)
                            %for each redundant value, get all the indices in objvals
                            redundantIdx = find(m.results(resultIdx).objvals(:,obj) == redundantVals(r));
                            if(length(redundantIdx) < 1), disp('Warning, something is wrong, expecting redundantIdx to be greater than 1'); end

                            %now since assuming rest of the objective values are already sorted, the last element in this list
                            %is the Pareto optimal, the rest would be weak Pareto optimal and are the ones to exclude
                            IdxToExclude = [IdxToExclude;redundantIdx(1:end-1)];
                        end
                        
                        %update m.results.objvals and m.results.sols
                        m.results(resultIdx).objvals = m.results(resultIdx).objvals(idx_all ~= IdxToExclude,:);
                        m.results(resultIdx).sols = m.results(resultIdx).sols(idx_all ~= IdxToExclude,:);
                    end
                end
            end
        end
        
    end
    
    methods (Access = public)
        %generate plots from m.solutions or given sol
        function generatePlotsPerSol(m, sol)
            if nargin < 2, sol = m.solutions; end
            solsN = size(sol.sols);
            regimens = 0:solsN(1)-1; %regimens = index starts with 0
            utopia_n = sol.utopia_n;
            nadir_n = sol.nadir_n;
            
            %% plot objective space
            figure;
            hold on;
            %plot all the points (if exist)
            if(~isempty(sol.all_cT))
                scatter(sol.all_cT(:,1),sol.all_cT(:,2), ...
                               'MarkerEdgeColor',[58/255 95/255 205/255]);
            end
            if(~isempty(sol.all_cF))
                scatter(sol.all_cF(:,1),sol.all_cF(:,2), ...
                               'MarkerEdgeColor',[112/255 128/255 144/255]);
            end
            
            %plot Pareto frontier points
            scatter(sol.objvals(:,1),sol.objvals(:,2), ...
                    'MarkerEdgeColor',[0/255 0/255 0/255], ...
                    'MarkerFaceColor',[220/255 28/255 28/255]);

            %plot compromise solution
            scatter(sol.objvals_compr(:,1),sol.objvals_compr(:,2), ...
                    'MarkerEdgeColor',[0/255 0/255 0/255], ...
                    'MarkerFaceColor',[20/255 140/255 20/255]);
            text(sol.objvals_compr(:,1)+0.02,sol.objvals_compr(:,2)+0.01,'Compromise');
            plot([sol.objvals_compr(:,1) utopia_n(1)],[sol.objvals_compr(:,2) utopia_n(2)],'--', ...
                    'Color',[112/255 128/255 144/255]);

            %plot utopia/nadir points
            scatter(utopia_n(1),utopia_n(2),'k','filled');
            text(utopia_n(1)+0.02,utopia_n(2)+0.05,'Utopia');
            scatter(nadir_n(1),nadir_n(2),'k','filled');
            text(nadir_n(1)+0.02,nadir_n(2)+0.05,'Nadir');

            %draw lines connecting nadir and utopia
            %position hardcoded here - assumes utopia=(1,0); nadir=(0,1)
            plot([0 1],[1 1],'--','Color',[112/255 128/255 144/255]);
            plot([1 1],[1 0],'--','Color',[112/255 128/255 144/255]);
            
            title('Objective space');
            xlabel(m.objNameMap(char(sol.objs(1))));
            ylabel(m.objNameMap(char(sol.objs(2))));
            %Utilities.applyPlotStyles(gcf);
            hold off;
            
            set(gca,'XLim',[0 1]);
            set(gca,'YLim',[0 1]);
            %m.saveFig(gcf);
            
            %% plot heatmap of all Pareto solutions
            sortsols = false; %sort based on drug frequency
            sortedIdx = 1:solsN(2);
            if(sortsols)
                freq = sum(sol.sols,1);
                [~,sortedIdx] = sort(freq,'descend');
            end
            figure;
            Utilities.heatmap_m(sol.sols(:,sortedIdx), ...
                            'collabels', m.drug_labels(sortedIdx), ...
                            'rowlabels', regimens, ...
                            'xlabel', 'Drugs', ...
                            'ylabel', 'Regimens', ...
                            'range', [0 1], ...
                            'rotateXLabel', 0, ...
                            'cmap', Utilities.getcmap('hot',true), ... %can also use gray_binary
                            'cbar', false);
            ht = title('Solution space');
            p = get(gcf, 'position');
            p(3) = p(3)*3;
            p(4) = p(4)*0.7;
            set(gcf,'position',p);
            Utilities.applyPlotStyles(gcf);
            %m.saveFig(gcf);

            %% plot histogram of drugs use
            xhistx = 1:solsN(2);
            xhist = sum(sol.sols,1);
            yhistx = regimens;
            yhist = sum(sol.sols,2)';
            
            figure;
            bar(xhistx,xhist,1);
            [m1 m2] = Utilities.barAxisRange(xhistx);
            set(gca,'XLim',[m1 m2]);
            xlabel('Drugs');
            ylabel('Number of regimens');
            set(gca,'XTick',1:m.drugsTotal, ...
                    'XTickLabel',m.drug_labels);
            %Utilities.applyPlotStyles(gcf);
            %m.saveFig(gcf);
            
            %% plot histogram of regimens use
            figure;
            bar(yhistx,yhist,1);
            [m1 m2] = Utilities.barAxisRange(yhistx);
            set(gca,'XLim',[m1 m2]);
            xlabel('Regimens');
            ylabel('Number of drugs');
            %Utilities.applyPlotStyles(gcf);
            %m.saveFig(gcf);

            %% plot heatmap of obj func values for all Pareto solutions
            figure;
            xn = 1;
            subplot(1,2,1);
            imagesc(sol.objvals(:,1), [min(sol.objvals(:,1)) max(sol.objvals(:,1))]);
            ylabel('Regimens');
            set(gca,'TickLength',[0 0], ...
                    'XTick',1:xn, ...
                    'XTickLabel',{m.objNameMap(char(sol.objs(1)))}, ...
                    'YTick',1:length(regimens), ...
                    'YTickLabel',regimens);
            colormap(hot);
            hcb = colorbar('location','eastoutside');
            set(hcb,'YTick',[0 0.5 1]);
            axis([0.5 xn+0.5 0.5 solsN(1)+0.5]);
            %Utilities.applyPlotStyles(gcf);
            
            subplot(1,2,2);
            imagesc(sol.objvals(:,2), [min(sol.objvals(:,2)) max(sol.objvals(:,2))]);
            ylabel('Regimens');
            set(gca,'TickLength',[0 0], ...
                    'XTick',1:xn, ...
                    'XTickLabel',{m.objNameMap(char(sol.objs(2)))}, ...
                    'YTick',1:length(regimens), ...
                    'YTickLabel',regimens);
            colormap(hot);
            hcb = colorbar('location','eastoutside');
            set(hcb,'YTick',[0 0.5 1]);
            axis([0.5 xn+0.5 0.5 solsN(1)+0.5]);
            %set(gca,'ydir','reverse'); %somehow Matlab reverses yaxis when using hold on on imagesc, flip it back
            %Utilities.applyPlotStyles(gcf);
            %m.saveFig(gcf);
            
            Utilities.applyPlotStylesToAll();
        end
        
        function generatePlots(m, anlyz)
            %generate all plots from m.analysis or given anlyz
            if nargin<2, anlyz = m.analysis; end
            
            if(isempty(anlyz))
                disp('Analysis data structure is empty. There is nothing to plot. Call m.analyze first.');
                return;
            end
            
            %% distribution of number of drugs
            figure;
            subplot(2,2,1);
            hold on;
            p = anlyz.drugsN/sum(anlyz.drugsN);
            b = bar(0:m.drugsTotal,p,1);
            set(b,'FaceColor',[100/255 100/255 100/255],'EdgeColor',[50/255 50/255 50/255]);
            errorbar(0:m.drugsTotal,p,sqrt(p.*(1-p)./length(m.results)),'.k','MarkerSize',2);
            set(gca,'XLim',[-0.5 m.drugsTotal+0.5]);
            xlabel('Number of drugs');
            ylabel('Count');
            title('Distribution of number of drugs');
            hold off;
            
            %% distribution of drug used
            subplot(2,2,2);
            hold on;
            p = anlyz.drugs/sum(anlyz.drugs);
            b = barh(1:m.drugsTotal,p,1);
            set(b,'FaceColor',[100/255 100/255 100/255],'EdgeColor',[50/255 50/255 50/255]);
            %h = herrorbar(p,1:m.drugsTotal,sqrt(p.*(1-p)./length(m.results)),'.k');
            %set(h,'MarkerSize',2);
            set(gca,'YLim',[0.5 m.drugsTotal+0.5]);
            ylabel('Drugs');
            xlabel('Count');
            title('Distribution of drug used');
            set(gca,'YTick',1:m.drugsTotal, ...
                    'YTickLabel',m.drug_labels);
            hold off;
            
            %% distribution of number of subpopulations simulated
            subplot(2,2,3);
            b = bar(1:m.subpopTotal,anlyz.subpopN/sum(anlyz.subpopN),1);
            set(b,'FaceColor',[100/255 100/255 100/255],'EdgeColor',[50/255 50/255 50/255]);
            set(gca,'XLim',[-0.5 m.subpopTotal+0.5]);
            xlabel('Number of subpopulations');
            ylabel('Count');
            title('Distribution of number of subpopulations');
            
            %% distribution of subpopulations used
            subplot(2,2,4);
            b = barh(1:m.subpopTotal,anlyz.subpop/sum(anlyz.subpop),1);
            set(b,'FaceColor',[100/255 100/255 100/255],'EdgeColor',[50/255 50/255 50/255]);
            set(gca,'YLim',[0.5 m.subpopTotal+0.5]);
            ylabel('Subpopulations');
            xlabel('Count');
            title('Distribution of subpopulations used');
            set(gca,'YTick',1:m.subpopTotal, ...
                    'YTickLabel',m.subpop_labels);
            %m.saveFig(gcf);
            
            Utilities.resizeFig(gcf,1.5,1.5);
            Utilities.applyPlotStyles(gcf);
            
            %% correlation between popN and drugN
            figure;
            scatter(anlyz.corr_popNdrugN(:,1),anlyz.corr_popNdrugN(:,2), 100);
            xlabel('Number of subpopulations');
            ylabel('Number of drugs');
            title('Correlation between popN and drugN');
            
            %% heatmap of drug usage vs number of subpopulations
            figure;
            Utilities.heatmap_m(anlyz.popNdrug_matrix, ...
                                'collabels', 1:m.subpopTotal, ...
                                'rowlabels', m.drug_labels, ...
                                'cmap', Utilities.getcmap('hot'), ...
                                'cbarlabel', 'Frequency', ...
                                'rotateXLabel', 0, ...
                                'xlabel', 'Number of subpopulations');
            title('Drug usage vs number of subpopulations');
            Utilities.resizeFig(gcf,1.5,1);
            Utilities.applyPlotStyles(gcf);
            
            figure;
            Utilities.heatmap_m(anlyz.corr_popNdrug', ...
                                'collabels', [], ...
                                'rowlabels', m.drug_labels, ...
                                'cmap', Utilities.getcmap('rwb'), ...
                                'cbarloc', 'eastoutside', ...
                                'cbarlabel', 'Spearman correlation', ...
                                'range', [-1 1], ...
                                'rotateXLabel', 0);
            Utilities.resizeFig(gcf,1/3.0,1);
            Utilities.applyPlotStyles(gcf);
            
            %% correlation of drug usage vs number of subpopulations
            figure;
            for i = 1:m.drugsTotal
                subplot(2,11,i);
                scatter(1:m.subpopTotal,anlyz.popNdrug_matrix(i,:));
                title({m.drug_labels{i},num2str(rounddec(anlyz.corr_popNdrug(i),10))});
                set(gca,'XLim',[0 m.subpopTotal]);
                
                if(sum(anlyz.popNdrug_matrix(i,:)) == 0)
                    set(gca,'YLim',[0 1]);
                end
            end
            p = get(gcf, 'position');
            Utilities.resizeFig(gcf,2.5,1/1.1);
            Utilities.applyPlotStyles(gcf,'small');
            
            %correlation freq and corr of # of drugs vs. # of subpopulations
            figure;
            scatter(sum(anlyz.corr_drug,1)',anlyz.corr_popNdrug', 100);
            c = corr(sum(anlyz.corr_drug,1)',anlyz.corr_popNdrug','type','Spearman','rows','complete'); %will ignore rows with NaNs
            title(['\rho = ' num2str(c)]);
            xlabel('Drug frequency');
            ylabel({'\rho (drug freq and number';'of subpopulations)'});
            Utilities.applyPlotStyles(gcf);
                            
            %% pairwise scatter plots - this is very time consuming
            %{
            figure;
            for i = 1:m.drugsTotal
                for j = 1:m.subpopTotal
                    subplot(m.drugsTotal, m.subpopTotal, (i-1)*m.subpopTotal+j);
                    scatter(m.analysis.corr_pop(:,j), m.analysis.corr_drug(:,i),20);
                    set(gca,'XTickLabel',[],'YTickLabel',[]);
                    if(i == 1), title(m.subpop_labels(j)); end
                    if(j == 1), ylabel(m.drug_labels(i)); end
                end
            end
            %}
            
            %% SENSITIVITY ANALYSIS
            % plot results with order of drugs based on frequency
            freq = anlyz.corr_variables(:,3);
            [freq_s freq_idx] = sort(freq);

            %% pairwise correlation plot
            % flip back the matrix, b/c matrix is flipped when using
            % imagesc (in heatmap_m)
            figure;
            Utilities.heatmap_m(flipud(anlyz.corr(freq_idx,:)), ...
                                'collabels', m.subpop_labels, ...
                                'rowlabels', flipud(anlyz.drugsIncluded_labels(freq_idx)), ...
                                'range', [-1 1], ...
                                'cbarlabel', 'Point-biserial correlation');
            Utilities.resizeFig(gcf,1.5,1);
            Utilities.applyPlotStyles(gcf);

            figure;
            Utilities.heatmap_m(flipud(abs(anlyz.corr(freq_idx,:))), ...
                                'collabels', m.subpop_labels, ...
                                'rowlabels', flipud(anlyz.drugsIncluded_labels(freq_idx)), ...
                                'range', [0 1], ...
                                'cmap', hot(100), ...
                                'cbarlabel', '|Point-biserial correlation|');
            Utilities.resizeFig(gcf,1.5,1);
            Utilities.applyPlotStyles(gcf);
            
            figure;
            boxplot(anlyz.corr(freq_idx,:)','orientation','horizontal', ...
                                'notch','off', ...
                                'symbol','ko', ...
                                'colors','k', ...
                                'labels',anlyz.drugsIncluded_labels(freq_idx));
            xlabel('Point-biserial correlation');
            Utilities.resizeFig(gcf,1,2);
            Utilities.applyPlotStyles(gcf);
            
            %% frequency of corr variables
            figure;
            for i = 1:5
                freq = anlyz.corr_variables(:,i);
                labels = m.drug_labels(anlyz.drugsToInclude);
                
                subplot(2,3,i);
                [freq_s idx] = sort(freq);
                barh(freq_s,1);
                set(gca,'YLim', [0.5 length(labels)+0.5], ...
                        'YTick', 1:length(labels), ...
                        'YTickLabel', labels(idx));
                xlabel('Value');
                ylabel('Drugs');
                title(anlyz.corr_variables_labels(i));
            end
            Utilities.resizeFig(gcf,2,1.5);
            Utilities.applyPlotStyles(gcf);
            
            %% correlation matrix            
            clustergram(anlyz.corr_matrix, ...
                        'RowLabels',anlyz.corr_variables_labels, ...
                        'ColumnLabels',anlyz.corr_variables_labels, ...
                        'colormap',Utilities.getcmap('rwb'), ...
                        'standardize',3, ...
                        'ShowDendrogram', 'on');

             %individual calls to Utilities.applyPlotStyles will already
             %also called Utilities.applyPlotBasicStyles. But apply
             %plotbasicstyles to all figure, just in case.
             Utilities.applyPlotBasicStylesToAll();
        end
        
        %% specific plot functions
        function plotDrugSubpop(m, drug, subpop, analyz)
            if nargin<4, anlyz = m.analysis; end
            
            if(isempty(anlyz))
                disp('Analysis data structure is empty. There is nothing to plot. Call m.analyze first.');
                return;
            end
            
            if(isempty(drug) || isempty(subpop))
                disp('Input drug or subpop is empty.');
                return;
            end
            
            figure;
            d_pre = find(strcmp([m.drug_labels],drug)); %index of drug, before excluding ones with freq = 0
            d = find(strcmp([anlyz.drugsIncluded_labels],drug)); %index of drug, after excluding ones with freq = 0
            p = find(strcmp([m.subpop_labels],subpop));
            if(~isempty(d) && ~isempty(p))
                scatter(anlyz.corr_pop(:,p), anlyz.corr_drug(:,d_pre), 100);
                set(gca,'YLim',[-0.5 1.5], 'YTick', [0 1]);
                title(['r_{pb} = ' num2str(anlyz.corr(d,p))]);
                ylabel(anlyz.drugsIncluded_labels(d));
                xlabel(m.subpop_labels(p));
            else
                disp(['Tried to plot sensitivity of drug vs. subpopulation, ', ...
                     'but can''t find specified drug/subpopulation.']);
            end
            Utilities.applyPlotStyles(gcf);
        end
    end
    
    methods (Access = private)
        function saveFig(m, h)
            m.figN = m.figN + 1;
            saveas(h,[m.outputsdir 'fig_' num2str(m.figN)],'epsc');
        end
        
        function sol = getSol(m, resultIdx, solIdx)
            if(isnumeric(solIdx))
                sol = m.results(resultIdx).sols(solIdx,:);
            else
                switch(solIdx)
                    case 'lowest'
                        sol = m.results(resultIdx).sols(1,:);
                    case 'highest'
                        sol = m.results(resultIdx).sols(end,:);
                    case 'compromise'
                        sol = m.results(resultIdx).sols_compr;
                end
            end
        end
        
    end
    
    %Analyses for given instances
    methods (Static)
        function [consensus, consensusSum, prop] = getConsensus(m0,m1,generatePlots)
            %Compares solutions between results m0 and m1 (with OR logic)
            %INPUT: 1) m0: simulate result to be compared to
            %           2) m1: vector of simulate results used to compare to m0
            %OUTPUT: 1) consensus data struct store drugs sols that are same/different beween m0 vs m1
            %                   NOTE: column index here starts with drugsN=1; results.sol starts with drugsN=0
            %                   exclude solIdx=1, no drug
            %                       row: pop struct
            %                       col: regimen (# of drugs in combo)
            %                       value: boolean, 0 disagree; 1 agree
            %               2) consensusSum = sum of consensus values
            %                       first row: # of pop that disagrees; second row: # of pop that agrees
            %                       col: regimen
            %               3) prop: proportions of popstructs that are different in drug combo that contains/not contain single-best drug (based on m1)
            %                       first row: # of pop that contains; second row: # of pop that does not contains
            %                       col: regimen
            
            if(nargin < 3), generatePlots = false; end
            consensus = -1*ones(length(m0.results),size(m0.results(1).sols,1)-1);

            %consensus matrix (popstruct versus regimen)
            for i = 1:length(m0.results)
                for j  = 2:size(m0.results(i).sols,1) %first solution is no-drug, exclude this
                    m0N = sum(m0.results(i).sols(j,:),2);
                    agree = false; %logic used for combining multiple m1: OR gate
                    for r = 1:length(m1)
                        m1N = sum(m1(r).results(i).sols(j,:),2);
                        if(m0N ~= m1N)
                            error('Number of drugs in combination is different for this instance between m0 and m1.');
                        end
                        agree = (agree | isequal(m0.results(i).sols(j,:),m1(r).results(i).sols(j,:)));
                    end
                    
                    %agree = isequal(m0.results(i).sols(j,:),m1.results(i).sols(j,:));
                    consensus(i, j-1) = agree;
                    
                end
            end

            %overall consensus over all populations
            consensusSum = -1*ones(2,size(consensus,2)); %row 1: disagree; row 2: agree; each col different drug combo
            consensusSum(2,:) = sum(consensus,1);
            consensusSum(1,:) = size(consensus,1)-consensusSum(2,:);
            
            %proportions of popstruct containing single-best drug for m1
            prop = zeros(2,size(consensus,2)); %row 1=contains single-best for predominant; row 2=does not contain...; col=regimens
            for regimen = 1:size(consensus,2) %loop through each regimen 
                for i = 1:size(consensus,1) %loop through each pop struct
                    if(~consensus(i,regimen)) %analyze only ones that drug sols differ
                        
                        contains = false; %logic used for combining multiple m1: OR gate
                        for r = 1:length(m1)
                            sol = m1(r).results(i).sols(2,:); %solution for 1-drug combination for m1
                            if(sum(sol) ~= 1), error('something is wrong, single-best is not a single drug'); end
                            drug_idx = find(sol == 1);
                            %note the consensus and sols index is offset by 1, since
                            %consensus excluded the drugsN=0 that's in sol
                            contains = (contains | (m0.results(i).sols(regimen+1,drug_idx) == 1));
                        end
                        
                        if(contains)
                            prop(1,regimen) = prop(1,regimen)+1;
                        else
                            prop(2,regimen) = prop(2,regimen)+1;
                        end
                    end
                end
            end
            
            %% Plots
            if(generatePlots)
                % plot stacked bar graph
                figure;
                h=bar(1:size(consensusSum,2),consensusSum',1,'stack');
                set(gca,'XLim',[0.5 size(consensusSum,2)+0.5], ...
                        'XTick',[1:size(consensusSum,2)], ...
                        'XTickLabel',1:size(consensusSum,2));
                legend('Different','Same');
                xlabel('Regimens');
                ylabel('Frequency');
                Utilities.applyPlotStyles(gcf);
                set(h(1),'Facecolor',[59/255 158/255 59/255]);
                set(h(2),'Facecolor',[50/255 108/255 191/255]);

                % compare what proportion of optimal drug combo contain also single-best drug for predominant
                for regimen = 1:size(prop,2) %loop through each regimen 
                    figure;
                    if(prop(1,regimen) == 0)
                        pie(prop(2,regimen));
                        %legend('does not contain');
                        colormap([50/255 50/255 50/255]);
                    elseif(prop(2,regimen) == 0)
                        pie(prop(1,regimen));
                        %legend('contains');
                        colormap([200/255 200/255 200/255]);
                    else
                        pie(prop(:,regimen));
                        %legend('contains','does not contain');
                        colormap([200/255 200/255 200/255;50/255 50/255 50/255]);
                    end
                end

                if(~(length(m1) >1))
                    % plot distribution of efficacy difference of the last combination
                    effDiff = [];
                    for i = 1:size(consensus,1)
                        %beware the index for consensus and sols is different; see description
                        %of consensus data struct above
                        if(~consensus(i,end)) %analyze only ones that drug sols differ
                            e1 = m0.results(i).sols(end,:)*m0.dEff*m0.results(i).pop;
                            e2 = m1.results(i).sols(end,:)*m0.dEff*m0.results(i).pop;
                            effDiff = [effDiff;e1-e2];
                        end
                    end
                    figure;
                    [f n] = hist(effDiff,30);
                    bar(n,f,1);
                    title('Last regimen');
                    ylabel('Frequency');
                    xlabel('Difference in efficacy');
                    Utilities.applyPlotStyles(gcf);
                end
            end
        end
        
        function results_combined = combineResults(results, resultsCombinedName)
            %combine multiple results
            %INPUT:	1) results: a vector of Modeling instances to combine
            %            2) resultsCombinedName: combined results name, to be set equal to m.resultsCombined
            %MODIFIES: this only modifies m.results; all rest of struct values are same as results(1)
            if nargin < 2, resultsCombinedName = 'combined'; end
            
            if(length(results) > 1)
                results_combined = results(1);
                for i = 2:length(results)
                    results_combined.results = [results_combined.results; results(i).results];
                end
                results_combined.resultsCombined = resultsCombinedName;
            else
                disp('Note: results only contain one data struct; no need to combine results.');
            end
        end
        
    end
end

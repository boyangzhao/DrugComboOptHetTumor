%% Utilities Class for combination therapy optimization on heterogeneous populations
% Author: Boyang Zhao
% Contains dependent static methods that other classes need

classdef Utilities < handle
    
    methods (Static)
        % CALCULATION METHODS
        function n = rounddec(n,d)
            n = round(n*10*d)/(10*d);
        end
        
        function t = comprval(x,y,tol)
            %Double precision floating point value comparsion
            if nargin < 3, tol = 1e-7; end
            t = abs(x-y)<tol;
        end
        
        function c = corr_pb(X, Y)
            %Calculates point-biserial correlation
            %Assume X is the binary data
            
            if(sum(isnan(X)) + sum(isnan(Y)) > 0)
                %return NaN is X or Y contains NaN
                c = nan;
            elseif(isempty(X) || isempty(Y))
                c = nan;
            else
                n0 = sum(X == 0);
                n1 = sum(X == 1);
                n = length(X);
                M0 = mean(Y(X == 0));
                M1 = mean(Y(X == 1));
                s = std(Y,1);

                c = ((M1-M0)/s)*sqrt((n0/n)*(n1/n));
            end
        end
        
        function mat = rotationMatrix(angle,degree)
            %angle in radians
            if nargin < 2, degree = false; end
            if degree, angle = degtorad(angle); end
            mat = [cos(angle), -1*sin(angle); sin(angle), cos(angle)];
        end
        
        function rcritical = calcRcritical(pval, df)
            %calculates the critical correlation coefficient value given
            %p-value threshold and degrees of freedom (df)
            tcritical = -1*tinv(1-pval/2,df);
            syms r;
            rcritical = eval(solve( 'r*sqrt(6-2)/sqrt(1-r^2) = tcritical',r));
            
            %check for errors
            temp = abs(rcritical);
            if(~comprval(temp(1),temp(2)))
                disp('Something went wrong with the rho critical calculation... solver solutions do not match');
                disp(['Solutions were ' num2str(rcritical(1)) ' and ' num2str(rcritical(1))]);
            end
            
            %the solutions should be the same, just opposite sign, so just need one of them
            rcritical = abs(rcritical(1));
        end
        
        % FIGURE GENERATION METHODS
        function [min max] = barAxisRange(x)
            if(length(x) == 1)
                %there's only one value in x, range will be +/- 0.5
                min = x - 0.5;
                max = x + 0.5;
            elseif(length(x) >= 2)
                w = (x(2)-x(1))/2;
                min = x(1)-w;
                max = x(end)+w;
            end
        end
        
        function resizeFig(h, width, height)
            %resize figure size (using position property)
            p = get(h, 'position');
            p(3) = p(3)*width;
            p(4) = p(4)*height;
            set(gcf,'position',p);
        end
        
        function applyPlotStylesToAll()
            h = get(0, 'Children');
            for i = 1:length(h)
                Utilities.applyPlotStyles(figure(h(i)));
            end
        end
        
        function applyPlotBasicStylesToAll()
            h = get(0, 'Children');
            for i = 1:length(h)
                Utilities.applyPlotBasicStyles(figure(h(i)));
            end
        end
        
        function applyPlotStyles(hf, plotSize)
            %input: hf = figure handle
            if(nargin < 2), plotSize = 'medium'; end
            switch(plotSize)
                case 'small'
                    markerSize = 20;
                    fontSize_label = 11;
                    fontSize_intext = 10;
                case 'medium'
                    markerSize = 100;
                    fontSize_label = 12;
                    fontSize_intext = 11;
                case 'large'
                    markerSize = 125;
                    fontSize_label = 13;
                    fontSize_intext = 12;
            end
            
            %apply styles independent of plot size
            Utilities.applyPlotBasicStyles(hf);
            
            %apply styles that may be plot size sensitive
            %e.g. font size, line width size, marker size, etc.
            ha = findobj(hf,'type','axes');
            for i = 1:length(ha)
                set(ha(i), 'FontSize', fontSize_label, ...
                           'LineWidth', 0.5);
                l = [get(ha(i),'XLabel') get(ha(i),'YLabel') get(ha(i),'ZLabel') get(ha(i),'Title')];    
                set(l, 'FontSize', fontSize_label);

                %group objects
                h2 = findobj(ha(i),'type','hggroup');
                try,set(h2, 'SizeData',markerSize);catch;end
                try,set(h2, 'LineWidth',1.5);catch;end
                try,set(h2, 'FaceColor',[100/255 100/255 100/255]);catch;end
                try,set(h2, 'EdgeColor',[50/255 50/255 50/255]);catch;end

                %line styles
                h2 = findobj(ha(i),'type','line');
                set(h2, 'LineWidth',0.5, ...
                        'Color','k');

                %text styles
                h2 = findobj(ha(i),'type','text');
                set(h2, 'FontSize',fontSize_intext, ...
                        'FontWeight','bold');
            end
        end
        
        function applyPlotBasicStyles(hf)
            %apply styles that are independent of the size of the plot
            %mainly: text font name = arial (normal weight), and color = black
            ha = findobj(hf,'type','axes');
            for i = 1:length(ha)
                set(ha(i), 'FontName', 'Arial', ...
                           'FontWeight', 'normal', ...
                           'Color', 'w', ...
                           'XColor', 'k', ...
                           'YColor', 'k', ...
                           'ZColor', 'k');
                l = [get(ha(i),'XLabel') get(ha(i),'YLabel') get(ha(i),'ZLabel') get(ha(i),'Title')];    
                set(l, 'FontName', 'Arial', ...
                       'FontWeight', 'normal', ...
                       'Color', 'k');

                %text styles
                h2 = findobj(ha(i),'type','text');
                set(h2, 'FontName', 'Arial', ...
                        'Color','k');
            end
        end
        
        function ht = rotateXLabel(h, angle)
            if nargin < 2, angle = 90; end

            %get locations
            txtlabels = get(h,'XTickLabel');
            xloc = get(h,'XTick');
            yloc = get(h,'XTick');

            if(~isempty(txtlabels))
                %erase current label and add in rotated label
                set(h,'XTickLabel',[]);
                ht = text(xloc, (yloc(1)-0.8*(yloc(2)-yloc(1)))*ones(1,size(txtlabels,1)), txtlabels, ...
                          'rotation', angle);

                %reposition figure to accommodate text
                p = get(h,'Position');
                p(2) = p(2)-p(4)*0.1;
                set(h,'Position',p);
            end
        end
        
        function heatmap_m(datamatrix, varargin)
            %Generates heatmap
            %dependencies: Utilities.getcmap, Utilities.rotateXLabel
            
            %defaults
            border = true;
            cmap = Utilities.getcmap('rwb');
            cbar = true;
            cbarloc = 'southoutside';
            cbarlabel = '';
            range = [min(min(datamatrix)) max(max(datamatrix))];
            rowlabels=[];
            collabels=[];
            rowlabelangle = 90;
            horzlabel = '';
            vertlabel = '';
            titlelabel = '';
            
            %get inputs
            n = size(varargin,2);
            i = 1;
            while(i <= n)
                switch(varargin{i})
                    case 'rowlabels'
                        rowlabels = varargin{i+1};
                    case 'collabels'
                        collabels = varargin{i+1};
                    case 'cbar'
                        cbar = varargin{i+1};
                    case 'cbarloc'
                        cbarloc = varargin{i+1};
                    case 'cbarlabel'
                        cbarlabel = varargin{i+1};
                    case 'cmap'
                        cmap = varargin{i+1};
                    case 'border'
                        border = varargin{i+1};
                    case 'range'
                        range = varargin{i+1};
                    case 'rotateXLabel'
                        rowlabelangle = varargin{i+1};
                    case 'xlabel'
                        horzlabel = varargin{i+1};
                    case 'ylabel'
                        vertlabel = varargin{i+1};
                    case 'title'
                        titlelabel = varargin{i+1};
                end
                i = i+2;
            end

            m = size(datamatrix,1);
            n = size(datamatrix,2);

            %plot
            imagesc(datamatrix, range);

            %styles/labels
            set(gca,'XTickLabel', collabels, ...
                        'XTick',1:n, ...
                        'YTickLabel', rowlabels, ...
                        'YTick',1:m, ...
                        'FontSize',11, ...
                        'TickLength',[0 0]);
            axis([0.5 n+0.5 0.5 m+0.5]);
            
            if(~isempty(horzlabel)), xlabel(horzlabel); end
            if(~isempty(vertlabel)), ylabel(vertlabel); end
            if(~isempty(titlelabel)), title(titlelabel); end
            
            if(~Utilities.comprval(rowlabelangle,0)), Utilities.rotateXLabel(gca, rowlabelangle); end
            
            %draw border
            if(border)
                hold on;
                for i = 0.5:n+0.5 %vertical lines
                    plot([i i], [0.5 m+0.5], '-', 'LineWidth', 1, 'color', [0 0 0]);
                end
                for i = 0.5:m+0.5 %horizontal lines
                    plot([0 n+0.5], [i i], '-', 'LineWidth', 1, 'color', [0 0 0]);
                end
                hold off;
            end

            %colormap/colorbar
            colormap(cmap);
            if cbar
                c = colorbar('location',cbarloc);
                if(~isempty(cbarlabel))
                    if(~isempty(strfind(cbarloc,'south')) || ~isempty(strfind(cbarloc,'north')))
                        set(get(c,'xlabel'),'string',cbarlabel);
                    elseif(~isempty(strfind(cbarloc,'east')) || ~isempty(strfind(cbarloc,'west')))
                        set(get(c,'ylabel'),'string',cbarlabel);
                    end
                end
            end
        end
        
        function c = getcmap(cmap,reverse,range)
            
            %Colors
            %wr: white -> red
            %rwb: blue -> white -> red
            %rbw2: rwb, lighter colors
            %hot: dark red -> yellow/white
            %gray_binary
            
            if(nargin < 2), reverse = false; end %reverse colormap
            if(nargin < 3), range = [0 1]; end %return a fraction of the colormap
            
            switch(cmap)
                case 'wr'
                    c = [0.9686    0.9686    0.9686
                        0.9696    0.9642    0.9611
                        0.9705    0.9598    0.9536
                        0.9715    0.9555    0.9460
                        0.9724    0.9511    0.9385
                        0.9733    0.9467    0.9310
                        0.9743    0.9423    0.9235
                        0.9752    0.9379    0.9159
                        0.9762    0.9335    0.9084
                        0.9771    0.9291    0.9009
                        0.9780    0.9247    0.8933
                        0.9790    0.9203    0.8858
                        0.9799    0.9159    0.8783
                        0.9809    0.9115    0.8707
                        0.9818    0.9071    0.8632
                        0.9827    0.9027    0.8557
                        0.9837    0.8984    0.8482
                        0.9846    0.8940    0.8406
                        0.9856    0.8896    0.8331
                        0.9865    0.8852    0.8256
                        0.9875    0.8808    0.8180
                        0.9884    0.8764    0.8105
                        0.9893    0.8720    0.8030
                        0.9903    0.8676    0.7955
                        0.9912    0.8632    0.7879
                        0.9922    0.8588    0.7804
                        0.9907    0.8504    0.7696
                        0.9893    0.8419    0.7587
                        0.9879    0.8334    0.7479
                        0.9865    0.8249    0.7371
                        0.9851    0.8165    0.7263
                        0.9837    0.8080    0.7155
                        0.9823    0.7995    0.7046
                        0.9809    0.7911    0.6938
                        0.9795    0.7826    0.6830
                        0.9780    0.7741    0.6722
                        0.9766    0.7656    0.6613
                        0.9752    0.7572    0.6505
                        0.9738    0.7487    0.6397
                        0.9724    0.7402    0.6289
                        0.9710    0.7318    0.6180
                        0.9696    0.7233    0.6072
                        0.9682    0.7148    0.5964
                        0.9667    0.7064    0.5856
                        0.9653    0.6979    0.5747
                        0.9639    0.6894    0.5639
                        0.9625    0.6809    0.5531
                        0.9611    0.6725    0.5423
                        0.9597    0.6640    0.5315
                        0.9583    0.6555    0.5206
                        0.9569    0.6471    0.5098
                        0.9522    0.6362    0.5015
                        0.9475    0.6254    0.4932
                        0.9427    0.6146    0.4849
                        0.9380    0.6038    0.4765
                        0.9333    0.5929    0.4682
                        0.9286    0.5821    0.4599
                        0.9239    0.5713    0.4516
                        0.9192    0.5605    0.4433
                        0.9145    0.5496    0.4350
                        0.9098    0.5388    0.4267
                        0.9051    0.5280    0.4184
                        0.9004    0.5172    0.4100
                        0.8957    0.5064    0.4017
                        0.8910    0.4955    0.3934
                        0.8863    0.4847    0.3851
                        0.8816    0.4739    0.3768
                        0.8769    0.4631    0.3685
                        0.8722    0.4522    0.3602
                        0.8675    0.4414    0.3518
                        0.8627    0.4306    0.3435
                        0.8580    0.4198    0.3352
                        0.8533    0.4089    0.3269
                        0.8486    0.3981    0.3186
                        0.8439    0.3873    0.3103
                        0.8392    0.3765    0.3020
                        0.8336    0.3652    0.2966
                        0.8279    0.3539    0.2913
                        0.8223    0.3426    0.2860
                        0.8166    0.3313    0.2806
                        0.8110    0.3200    0.2753
                        0.8053    0.3087    0.2700
                        0.7997    0.2974    0.2646
                        0.7940    0.2861    0.2593
                        0.7884    0.2748    0.2540
                        0.7827    0.2635    0.2486
                        0.7771    0.2522    0.2433
                        0.7715    0.2409    0.2380
                        0.7658    0.2296    0.2326
                        0.7602    0.2184    0.2273
                        0.7545    0.2071    0.2220
                        0.7489    0.1958    0.2166
                        0.7432    0.1845    0.2113
                        0.7376    0.1732    0.2060
                        0.7319    0.1619    0.2006
                        0.7263    0.1506    0.1953
                        0.7206    0.1393    0.1900
                        0.7150    0.1280    0.1846
                        0.7093    0.1167    0.1793
                        0.7037    0.1054    0.1740
                        0.6980    0.0941    0.1686
                        0.6863    0.0904    0.1667
                        0.6745    0.0866    0.1649
                        0.6627    0.0828    0.1630
                        0.6510    0.0791    0.1611
                        0.6392    0.0753    0.1592
                        0.6275    0.0715    0.1573
                        0.6157    0.0678    0.1555
                        0.6039    0.0640    0.1536
                        0.5922    0.0602    0.1517
                        0.5804    0.0565    0.1498
                        0.5686    0.0527    0.1479
                        0.5569    0.0489    0.1460
                        0.5451    0.0452    0.1442
                        0.5333    0.0414    0.1423
                        0.5216    0.0376    0.1404
                        0.5098    0.0339    0.1385
                        0.4980    0.0301    0.1366
                        0.4863    0.0264    0.1347
                        0.4745    0.0226    0.1329
                        0.4627    0.0188    0.1310
                        0.4510    0.0151    0.1291
                        0.4392    0.0113    0.1272
                        0.4275    0.0075    0.1253
                        0.4157    0.0038    0.1235
                        0.4039         0    0.1216];
                case 'rwb'
                    c = [0.0196    0.1882    0.3804
                        0.0240    0.1967    0.3922
                        0.0284    0.2052    0.4039
                        0.0328    0.2136    0.4157
                        0.0372    0.2221    0.4275
                        0.0416    0.2306    0.4392
                        0.0460    0.2391    0.4510
                        0.0504    0.2475    0.4627
                        0.0547    0.2560    0.4745
                        0.0591    0.2645    0.4863
                        0.0635    0.2729    0.4980
                        0.0679    0.2814    0.5098
                        0.0723    0.2899    0.5216
                        0.0767    0.2984    0.5333
                        0.0811    0.3068    0.5451
                        0.0855    0.3153    0.5569
                        0.0899    0.3238    0.5686
                        0.0943    0.3322    0.5804
                        0.0987    0.3407    0.5922
                        0.1031    0.3492    0.6039
                        0.1075    0.3576    0.6157
                        0.1118    0.3661    0.6275
                        0.1162    0.3746    0.6392
                        0.1206    0.3831    0.6510
                        0.1250    0.3915    0.6627
                        0.1294    0.4000    0.6745
                        0.1347    0.4071    0.6781
                        0.1401    0.4141    0.6817
                        0.1454    0.4212    0.6853
                        0.1507    0.4282    0.6889
                        0.1561    0.4353    0.6925
                        0.1614    0.4424    0.6962
                        0.1667    0.4494    0.6998
                        0.1721    0.4565    0.7034
                        0.1774    0.4635    0.7070
                        0.1827    0.4706    0.7106
                        0.1881    0.4776    0.7142
                        0.1934    0.4847    0.7178
                        0.1987    0.4918    0.7214
                        0.2041    0.4988    0.7250
                        0.2094    0.5059    0.7286
                        0.2147    0.5129    0.7322
                        0.2201    0.5200    0.7358
                        0.2254    0.5271    0.7395
                        0.2307    0.5341    0.7431
                        0.2361    0.5412    0.7467
                        0.2414    0.5482    0.7503
                        0.2467    0.5553    0.7539
                        0.2521    0.5624    0.7575
                        0.2574    0.5694    0.7611
                        0.2627    0.5765    0.7647
                        0.2751    0.5843    0.7689
                        0.2875    0.5922    0.7732
                        0.2999    0.6000    0.7774
                        0.3123    0.6078    0.7816
                        0.3247    0.6157    0.7859
                        0.3371    0.6235    0.7901
                        0.3495    0.6314    0.7944
                        0.3619    0.6392    0.7986
                        0.3743    0.6471    0.8028
                        0.3867    0.6549    0.8071
                        0.3991    0.6627    0.8113
                        0.4115    0.6706    0.8155
                        0.4238    0.6784    0.8198
                        0.4362    0.6863    0.8240
                        0.4486    0.6941    0.8282
                        0.4610    0.7020    0.8325
                        0.4734    0.7098    0.8367
                        0.4858    0.7176    0.8409
                        0.4982    0.7255    0.8452
                        0.5106    0.7333    0.8494
                        0.5230    0.7412    0.8536
                        0.5354    0.7490    0.8579
                        0.5478    0.7569    0.8621
                        0.5602    0.7647    0.8664
                        0.5725    0.7725    0.8706
                        0.5824    0.7776    0.8734
                        0.5923    0.7826    0.8762
                        0.6022    0.7876    0.8791
                        0.6121    0.7926    0.8819
                        0.6220    0.7976    0.8847
                        0.6318    0.8027    0.8875
                        0.6417    0.8077    0.8904
                        0.6516    0.8127    0.8932
                        0.6615    0.8177    0.8960
                        0.6714    0.8227    0.8988
                        0.6813    0.8278    0.9016
                        0.6911    0.8328    0.9045
                        0.7010    0.8378    0.9073
                        0.7109    0.8428    0.9101
                        0.7208    0.8478    0.9129
                        0.7307    0.8529    0.9158
                        0.7405    0.8579    0.9186
                        0.7504    0.8629    0.9214
                        0.7603    0.8679    0.9242
                        0.7702    0.8729    0.9271
                        0.7801    0.8780    0.9299
                        0.7900    0.8830    0.9327
                        0.7998    0.8880    0.9355
                        0.8097    0.8930    0.9384
                        0.8196    0.8980    0.9412
                        0.8256    0.9009    0.9423
                        0.8315    0.9037    0.9434
                        0.8375    0.9065    0.9445
                        0.8435    0.9093    0.9456
                        0.8494    0.9122    0.9467
                        0.8554    0.9150    0.9478
                        0.8613    0.9178    0.9489
                        0.8673    0.9206    0.9500
                        0.8733    0.9235    0.9511
                        0.8792    0.9263    0.9522
                        0.8852    0.9291    0.9533
                        0.8911    0.9319    0.9544
                        0.8971    0.9347    0.9555
                        0.9031    0.9376    0.9565
                        0.9090    0.9404    0.9576
                        0.9150    0.9432    0.9587
                        0.9209    0.9460    0.9598
                        0.9269    0.9489    0.9609
                        0.9329    0.9517    0.9620
                        0.9388    0.9545    0.9631
                        0.9448    0.9573    0.9642
                        0.9507    0.9602    0.9653
                        0.9567    0.9630    0.9664
                        0.9627    0.9658    0.9675
                        0.9686    0.9686    0.9686
                        0.9696    0.9642    0.9611
                        0.9705    0.9598    0.9536
                        0.9715    0.9555    0.9460
                        0.9724    0.9511    0.9385
                        0.9733    0.9467    0.9310
                        0.9743    0.9423    0.9235
                        0.9752    0.9379    0.9159
                        0.9762    0.9335    0.9084
                        0.9771    0.9291    0.9009
                        0.9780    0.9247    0.8933
                        0.9790    0.9203    0.8858
                        0.9799    0.9159    0.8783
                        0.9809    0.9115    0.8707
                        0.9818    0.9071    0.8632
                        0.9827    0.9027    0.8557
                        0.9837    0.8984    0.8482
                        0.9846    0.8940    0.8406
                        0.9856    0.8896    0.8331
                        0.9865    0.8852    0.8256
                        0.9875    0.8808    0.8180
                        0.9884    0.8764    0.8105
                        0.9893    0.8720    0.8030
                        0.9903    0.8676    0.7955
                        0.9912    0.8632    0.7879
                        0.9922    0.8588    0.7804
                        0.9907    0.8504    0.7696
                        0.9893    0.8419    0.7587
                        0.9879    0.8334    0.7479
                        0.9865    0.8249    0.7371
                        0.9851    0.8165    0.7263
                        0.9837    0.8080    0.7155
                        0.9823    0.7995    0.7046
                        0.9809    0.7911    0.6938
                        0.9795    0.7826    0.6830
                        0.9780    0.7741    0.6722
                        0.9766    0.7656    0.6613
                        0.9752    0.7572    0.6505
                        0.9738    0.7487    0.6397
                        0.9724    0.7402    0.6289
                        0.9710    0.7318    0.6180
                        0.9696    0.7233    0.6072
                        0.9682    0.7148    0.5964
                        0.9667    0.7064    0.5856
                        0.9653    0.6979    0.5747
                        0.9639    0.6894    0.5639
                        0.9625    0.6809    0.5531
                        0.9611    0.6725    0.5423
                        0.9597    0.6640    0.5315
                        0.9583    0.6555    0.5206
                        0.9569    0.6471    0.5098
                        0.9522    0.6362    0.5015
                        0.9475    0.6254    0.4932
                        0.9427    0.6146    0.4849
                        0.9380    0.6038    0.4765
                        0.9333    0.5929    0.4682
                        0.9286    0.5821    0.4599
                        0.9239    0.5713    0.4516
                        0.9192    0.5605    0.4433
                        0.9145    0.5496    0.4350
                        0.9098    0.5388    0.4267
                        0.9051    0.5280    0.4184
                        0.9004    0.5172    0.4100
                        0.8957    0.5064    0.4017
                        0.8910    0.4955    0.3934
                        0.8863    0.4847    0.3851
                        0.8816    0.4739    0.3768
                        0.8769    0.4631    0.3685
                        0.8722    0.4522    0.3602
                        0.8675    0.4414    0.3518
                        0.8627    0.4306    0.3435
                        0.8580    0.4198    0.3352
                        0.8533    0.4089    0.3269
                        0.8486    0.3981    0.3186
                        0.8439    0.3873    0.3103
                        0.8392    0.3765    0.3020
                        0.8336    0.3652    0.2966
                        0.8279    0.3539    0.2913
                        0.8223    0.3426    0.2860
                        0.8166    0.3313    0.2806
                        0.8110    0.3200    0.2753
                        0.8053    0.3087    0.2700
                        0.7997    0.2974    0.2646
                        0.7940    0.2861    0.2593
                        0.7884    0.2748    0.2540
                        0.7827    0.2635    0.2486
                        0.7771    0.2522    0.2433
                        0.7715    0.2409    0.2380
                        0.7658    0.2296    0.2326
                        0.7602    0.2184    0.2273
                        0.7545    0.2071    0.2220
                        0.7489    0.1958    0.2166
                        0.7432    0.1845    0.2113
                        0.7376    0.1732    0.2060
                        0.7319    0.1619    0.2006
                        0.7263    0.1506    0.1953
                        0.7206    0.1393    0.1900
                        0.7150    0.1280    0.1846
                        0.7093    0.1167    0.1793
                        0.7037    0.1054    0.1740
                        0.6980    0.0941    0.1686
                        0.6863    0.0904    0.1667
                        0.6745    0.0866    0.1649
                        0.6627    0.0828    0.1630
                        0.6510    0.0791    0.1611
                        0.6392    0.0753    0.1592
                        0.6275    0.0715    0.1573
                        0.6157    0.0678    0.1555
                        0.6039    0.0640    0.1536
                        0.5922    0.0602    0.1517
                        0.5804    0.0565    0.1498
                        0.5686    0.0527    0.1479
                        0.5569    0.0489    0.1460
                        0.5451    0.0452    0.1442
                        0.5333    0.0414    0.1423
                        0.5216    0.0376    0.1404
                        0.5098    0.0339    0.1385
                        0.4980    0.0301    0.1366
                        0.4863    0.0264    0.1347
                        0.4745    0.0226    0.1329
                        0.4627    0.0188    0.1310
                        0.4510    0.0151    0.1291
                        0.4392    0.0113    0.1272
                        0.4275    0.0075    0.1253
                        0.4157    0.0038    0.1235
                        0.4039         0    0.1216];
                case 'hot'
                    c = [0.0104         0         0
                        0.0208         0         0
                        0.0313         0         0
                        0.0417         0         0
                        0.0521         0         0
                        0.0625         0         0
                        0.0729         0         0
                        0.0833         0         0
                        0.0938         0         0
                        0.1042         0         0
                        0.1146         0         0
                        0.1250         0         0
                        0.1354         0         0
                        0.1458         0         0
                        0.1563         0         0
                        0.1667         0         0
                        0.1771         0         0
                        0.1875         0         0
                        0.1979         0         0
                        0.2083         0         0
                        0.2188         0         0
                        0.2292         0         0
                        0.2396         0         0
                        0.2500         0         0
                        0.2604         0         0
                        0.2708         0         0
                        0.2813         0         0
                        0.2917         0         0
                        0.3021         0         0
                        0.3125         0         0
                        0.3229         0         0
                        0.3333         0         0
                        0.3438         0         0
                        0.3542         0         0
                        0.3646         0         0
                        0.3750         0         0
                        0.3854         0         0
                        0.3958         0         0
                        0.4063         0         0
                        0.4167         0         0
                        0.4271         0         0
                        0.4375         0         0
                        0.4479         0         0
                        0.4583         0         0
                        0.4688         0         0
                        0.4792         0         0
                        0.4896         0         0
                        0.5000         0         0
                        0.5104         0         0
                        0.5208         0         0
                        0.5313         0         0
                        0.5417         0         0
                        0.5521         0         0
                        0.5625         0         0
                        0.5729         0         0
                        0.5833         0         0
                        0.5938         0         0
                        0.6042         0         0
                        0.6146         0         0
                        0.6250         0         0
                        0.6354         0         0
                        0.6458         0         0
                        0.6563         0         0
                        0.6667         0         0
                        0.6771         0         0
                        0.6875         0         0
                        0.6979         0         0
                        0.7083         0         0
                        0.7188         0         0
                        0.7292         0         0
                        0.7396         0         0
                        0.7500         0         0
                        0.7604         0         0
                        0.7708         0         0
                        0.7813         0         0
                        0.7917         0         0
                        0.8021         0         0
                        0.8125         0         0
                        0.8229         0         0
                        0.8333         0         0
                        0.8438         0         0
                        0.8542         0         0
                        0.8646         0         0
                        0.8750         0         0
                        0.8854         0         0
                        0.8958         0         0
                        0.9063         0         0
                        0.9167         0         0
                        0.9271         0         0
                        0.9375         0         0
                        0.9479         0         0
                        0.9583         0         0
                        0.9688         0         0
                        0.9792         0         0
                        0.9896         0         0
                        1.0000         0         0
                        1.0000    0.0104         0
                        1.0000    0.0208         0
                        1.0000    0.0313         0
                        1.0000    0.0417         0
                        1.0000    0.0521         0
                        1.0000    0.0625         0
                        1.0000    0.0729         0
                        1.0000    0.0833         0
                        1.0000    0.0938         0
                        1.0000    0.1042         0
                        1.0000    0.1146         0
                        1.0000    0.1250         0
                        1.0000    0.1354         0
                        1.0000    0.1458         0
                        1.0000    0.1563         0
                        1.0000    0.1667         0
                        1.0000    0.1771         0
                        1.0000    0.1875         0
                        1.0000    0.1979         0
                        1.0000    0.2083         0
                        1.0000    0.2188         0
                        1.0000    0.2292         0
                        1.0000    0.2396         0
                        1.0000    0.2500         0
                        1.0000    0.2604         0
                        1.0000    0.2708         0
                        1.0000    0.2813         0
                        1.0000    0.2917         0
                        1.0000    0.3021         0
                        1.0000    0.3125         0
                        1.0000    0.3229         0
                        1.0000    0.3333         0
                        1.0000    0.3438         0
                        1.0000    0.3542         0
                        1.0000    0.3646         0
                        1.0000    0.3750         0
                        1.0000    0.3854         0
                        1.0000    0.3958         0
                        1.0000    0.4063         0
                        1.0000    0.4167         0
                        1.0000    0.4271         0
                        1.0000    0.4375         0
                        1.0000    0.4479         0
                        1.0000    0.4583         0
                        1.0000    0.4688         0
                        1.0000    0.4792         0
                        1.0000    0.4896         0
                        1.0000    0.5000         0
                        1.0000    0.5104         0
                        1.0000    0.5208         0
                        1.0000    0.5313         0
                        1.0000    0.5417         0
                        1.0000    0.5521         0
                        1.0000    0.5625         0
                        1.0000    0.5729         0
                        1.0000    0.5833         0
                        1.0000    0.5938         0
                        1.0000    0.6042         0
                        1.0000    0.6146         0
                        1.0000    0.6250         0
                        1.0000    0.6354         0
                        1.0000    0.6458         0
                        1.0000    0.6563         0
                        1.0000    0.6667         0
                        1.0000    0.6771         0
                        1.0000    0.6875         0
                        1.0000    0.6979         0
                        1.0000    0.7083         0
                        1.0000    0.7188         0
                        1.0000    0.7292         0
                        1.0000    0.7396         0
                        1.0000    0.7500         0
                        1.0000    0.7604         0
                        1.0000    0.7708         0
                        1.0000    0.7813         0
                        1.0000    0.7917         0
                        1.0000    0.8021         0
                        1.0000    0.8125         0
                        1.0000    0.8229         0
                        1.0000    0.8333         0
                        1.0000    0.8438         0
                        1.0000    0.8542         0
                        1.0000    0.8646         0
                        1.0000    0.8750         0
                        1.0000    0.8854         0
                        1.0000    0.8958         0
                        1.0000    0.9063         0
                        1.0000    0.9167         0
                        1.0000    0.9271         0
                        1.0000    0.9375         0
                        1.0000    0.9479         0
                        1.0000    0.9583         0
                        1.0000    0.9688         0
                        1.0000    0.9792         0
                        1.0000    0.9896         0
                        1.0000    1.0000         0
                        1.0000    1.0000    0.0156
                        1.0000    1.0000    0.0313
                        1.0000    1.0000    0.0469
                        1.0000    1.0000    0.0625
                        1.0000    1.0000    0.0781
                        1.0000    1.0000    0.0938
                        1.0000    1.0000    0.1094
                        1.0000    1.0000    0.1250
                        1.0000    1.0000    0.1406
                        1.0000    1.0000    0.1563
                        1.0000    1.0000    0.1719
                        1.0000    1.0000    0.1875
                        1.0000    1.0000    0.2031
                        1.0000    1.0000    0.2188
                        1.0000    1.0000    0.2344
                        1.0000    1.0000    0.2500
                        1.0000    1.0000    0.2656
                        1.0000    1.0000    0.2813
                        1.0000    1.0000    0.2969
                        1.0000    1.0000    0.3125
                        1.0000    1.0000    0.3281
                        1.0000    1.0000    0.3438
                        1.0000    1.0000    0.3594
                        1.0000    1.0000    0.3750
                        1.0000    1.0000    0.3906
                        1.0000    1.0000    0.4063
                        1.0000    1.0000    0.4219
                        1.0000    1.0000    0.4375
                        1.0000    1.0000    0.4531
                        1.0000    1.0000    0.4688
                        1.0000    1.0000    0.4844
                        1.0000    1.0000    0.5000
                        1.0000    1.0000    0.5156
                        1.0000    1.0000    0.5313
                        1.0000    1.0000    0.5469
                        1.0000    1.0000    0.5625
                        1.0000    1.0000    0.5781
                        1.0000    1.0000    0.5938
                        1.0000    1.0000    0.6094
                        1.0000    1.0000    0.6250
                        1.0000    1.0000    0.6406
                        1.0000    1.0000    0.6563
                        1.0000    1.0000    0.6719
                        1.0000    1.0000    0.6875
                        1.0000    1.0000    0.7031
                        1.0000    1.0000    0.7188
                        1.0000    1.0000    0.7344
                        1.0000    1.0000    0.7500
                        1.0000    1.0000    0.7656
                        1.0000    1.0000    0.7813
                        1.0000    1.0000    0.7969
                        1.0000    1.0000    0.8125
                        1.0000    1.0000    0.8281
                        1.0000    1.0000    0.8438
                        1.0000    1.0000    0.8594
                        1.0000    1.0000    0.8750
                        1.0000    1.0000    0.8906
                        1.0000    1.0000    0.9063
                        1.0000    1.0000    0.9219
                        1.0000    1.0000    0.9375
                        1.0000    1.0000    0.9531
                        1.0000    1.0000    0.9688
                        1.0000    1.0000    0.9844
                        1.0000    1.0000    1.0000];
                case 'rwb2'
                    c = [0.0000	0.0000	1.0000
                        0.0323	0.0323	1.0000
                        0.0645	0.0645	1.0000
                        0.0968	0.0968	1.0000
                        0.1290	0.1290	1.0000
                        0.1613	0.1613	1.0000
                        0.1935	0.1935	1.0000
                        0.2258	0.2258	1.0000
                        0.2581	0.2581	1.0000
                        0.2903	0.2903	1.0000
                        0.3226	0.3226	1.0000
                        0.3548	0.3548	1.0000
                        0.3871	0.3871	1.0000
                        0.4194	0.4194	1.0000
                        0.4516	0.4516	1.0000
                        0.4839	0.4839	1.0000
                        0.5161	0.5161	1.0000
                        0.5484	0.5484	1.0000
                        0.5806	0.5806	1.0000
                        0.6129	0.6129	1.0000
                        0.6452	0.6452	1.0000
                        0.6774	0.6774	1.0000
                        0.7097	0.7097	1.0000
                        0.7419	0.7419	1.0000
                        0.7742	0.7742	1.0000
                        0.8065	0.8065	1.0000
                        0.8387	0.8387	1.0000
                        0.8710	0.8710	1.0000
                        0.9032	0.9032	1.0000
                        0.9355	0.9355	1.0000
                        0.9677	0.9677	1.0000
                        1.0000	1.0000	1.0000
                        1.0000	1.0000	1.0000
                        1.0000	0.9677	0.9677
                        1.0000	0.9355	0.9355
                        1.0000	0.9032	0.9032
                        1.0000	0.8710	0.8710
                        1.0000	0.8387	0.8387
                        1.0000	0.8065	0.8065
                        1.0000	0.7742	0.7742
                        1.0000	0.7419	0.7419
                        1.0000	0.7097	0.7097
                        1.0000	0.6774	0.6774
                        1.0000	0.6452	0.6452
                        1.0000	0.6129	0.6129
                        1.0000	0.5806	0.5806
                        1.0000	0.5484	0.5484
                        1.0000	0.5161	0.5161
                        1.0000	0.4839	0.4839
                        1.0000	0.4516	0.4516
                        1.0000	0.4194	0.4194
                        1.0000	0.3871	0.3871
                        1.0000	0.3548	0.3548
                        1.0000	0.3226	0.3226
                        1.0000	0.2903	0.2903
                        1.0000	0.2581	0.2581
                        1.0000	0.2258	0.2258
                        1.0000	0.1935	0.1935
                        1.0000	0.1613	0.1613
                        1.0000	0.1290	0.1290
                        1.0000	0.0968	0.0968
                        1.0000	0.0645	0.0645
                        1.0000	0.0323	0.0323
                        1.0000	0.0000	0.0000];
                case 'gray_binary'
                    c = [1 1 1; 100/255 100/255 100/255];
            end
            
            if(reverse), c = flipud(c); end
            if ~isequal(range, [0,1])
                select = ceil(size(c,1)*range(1))+1:ceil(size(c,1)*range(2));
                c = c(select,:);
            end
        end
        
    end
end


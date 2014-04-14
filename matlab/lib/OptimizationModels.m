%% Optimization Models Class for combination therapy optimization on heterogeneous populations
% Author: Boyang Zhao
% Requires: CPlex libraries for Matlab

classdef OptimizationModels < handle & Modeling
    %settings
    properties (SetAccess = public, GetAccess = public)
        optModel;
    end
    
    %derived variables
    properties (SetAccess = protected, GetAccess = public)
        optModelName;
    end
    
    %results
    properties (SetAccess = private, GetAccess = protected)
        solutions;
    end
    
    %internal variables
    properties (SetAccess = private, GetAccess = protected)
        prevpts, optModelMap, objNameMap;
    end
    
    methods (Access = public)
        function m = OptimizationModels(dsets)
            %Constructor
            m = m@Modeling(dsets);
            m.prevpts = [];
            m.solutions = struct('objs', [], ...
                                 'pop', [], ...
                                 'utopia', [], ...
                                 'nadir', [], ...
                                 'utopia_n', [], ...
                                 'nadir_n', [], ...
                                 'normfactor', [], ...
                                 'sols', [], ...
                                 'objvals', [], ...
                                 'sols_compr', [], ... %compromise sol
                                 'objvals_compr', [], ... %obj values for compromise sol
                                 'all_cT', [], ... %all_cT: satifies constraint (cT = constraint true)
                                 'all_cF', []); %all_cF: does not satify constraint (cF = constraint false)
            m.optModel = '';
            objModelLabel_short = {'efftox'};
            objModelLabel_long = {'efftox_awT'};
            m.optModelMap = containers.Map(objModelLabel_short,objModelLabel_long);
            
            objLabel_short = {'eff','tox'};
            objLabel_long = {'Efficacy','Toxicity'};
            m.objNameMap = containers.Map(objLabel_short,objLabel_long);
            
            %check if Cplex library exist
            if ~exist('Cplex','class')
                %if not, try to add it to path if there's a cplex folder in
                %the current working directory
                if exist('cplex','dir')
                    addpath('cplex');
                else
                    error('Cannot find Cplex library. Unable to run OptimizationModels without Cplex');
                end
            end
        end
        
        function optimize(m)
            %Public optimization wrapper function
            %MODIFIES: m.solutions
            
            m.solutions.pop = m.popstruct_e;
            
            if m.changeToSubpopOnly > 0
                %only optimize based on the nth largest subpopulation
                %change m.popstruct_e so now it only contains 100% of
                %the nth largest subpopulation
                [val, idx] = sort(m.popstruct_e,'descend');
                
                %replace the population structure with 100% of the nth largest subpopulation
                idx_toSetFull = idx(m.changeToSubpopOnly);
                m.popstruct_e(idx_toSetFull) = 1;
                m.popstruct_e(setdiff(idx,idx_toSetFull)) = 0;
            end
            
            switch m.optModel
                case 'efftox' %efficacy + frontline + toxicity
                    m.solutions.objs = {'eff';'tox'};
                    if(m.listAllSolutions), m.findAllSolutions(); end
                    m.awTchebycheffOptWrapper(@m.efftox);
            end
            
            %clear the tracking variable (on if utopia/nadir has previously
            %been calculated
            m.prevpts = [];
        end
    end
    
    methods (Access = private)
%% Wrapper Functions for Optimization
        function awTchebycheffOptWrapper(m, optfunc)
            %wrapper function for augmented weighted Tchebycheff optimization
            rho = 10e-6;
            weights = 0:0.01:1;
            weights = [weights' (1-weights')];
            weights = [weights; 1 1]; %for compromise solution

            m.findPts();
            sols = -1*zeros(size(weights,1),m.drugsTotal); %preallocate
            for i = 1:size(weights,1)
                d = optfunc(rho,weights(i,:));
                sols(i,:) = d';
            end

            %only save the unique solutions
            m.solutions.sols = unique(sols(1:end,:),'rows');
            m.solutions.objvals = [m.calcObjVal(m.solutions.sols',m.solutions.objs(1)) ...
                                            m.calcObjVal(m.solutions.sols',m.solutions.objs(2))];
            
            %sort according to the objective 1 values
            [~, idx] = sort(m.solutions.objvals(:,1));
            m.solutions.objvals = m.solutions.objvals(idx,:);
            m.solutions.sols = m.solutions.sols(idx,:);

            %compromise solution
            m.solutions.sols_compr = sols(end,:);
            m.solutions.objvals_compr = [m.calcObjVal(m.solutions.sols_compr',m.solutions.objs(1)) ...
                                                    m.calcObjVal(m.solutions.sols_compr',m.solutions.objs(2))];
        end
        
        function findPts(m)
            %find utopia and nadir points
            
            if(length(m.solutions.objs) > 2)
                error('A maximum of two objectives allowed for determining utopia and nadir points');
            end
            
            %if already found the points for objectives obj, then don't
            %have to run the optimization again
            if(~isempty(m.prevpts) && isequal(m.prevpts,m.solutions.objs)), return; end
            
            m.solutions.utopia = [0 0];
            m.solutions.nadir = [0 0];
            m.solutions.utopia_n = m.solutions.utopia;
            m.solutions.nadir_n = m.solutions.utopia;
            
            for n = 1:length(m.solutions.objs)
                switch(char(m.solutions.objs(n)))
                    case 'eff' %objective 1 efficacy
                        m.eff(n);
                    
                    case 'tox' %objective 2 toxicity
                        m.tox(n);
                end
            end
            
            m.solutions.normfactor = 1./(m.solutions.utopia - m.solutions.nadir);
            m.prevpts = m.solutions.objs;
        end
        
%% Objective Function Calculations
        function val = calcObjVal(m, x, obj)
            if(isinteger(x))
                x = double(x);
            end
            switch(char(obj))
                case 'eff' %efficacy
                    val = m.normalize(x'*m.dEff*m.popstruct_e,1);
                case 'tox' %toxicity
                    val = m.normalize(-x'*m.dTox*ones(size(m.dTox,2),1),2,true);
            end
        end
        
        function x = normalize(m, x, idx, flip)
            %assumes maximization problem, with utopia >= x and nadir <= x
            %normalizes data to range 0 to 1
            %if flip is true, then flip the values (used for maximizing negative of obj and need to flip it back)
            if nargin <= 3, flip = false; end
            x = (x - m.solutions.nadir(idx))/(m.solutions.utopia(idx) - m.solutions.nadir(idx));
            if flip, x=-x+1; end
        end
        
        function findAllSolutions(m)
            %Note: can preallocate allv_cT and allv_cF
            m.findPts();

            if(Utilities.comprval(m.drugsN,-1))
                drugNRange = 0:m.drugsTotal;
            else
                drugNRange = m.drugsN;
            end
            for j = drugNRange
                if Utilities.comprval(j,0)
                    allx = -1; %mark as a non-empty vector, so this can still go through the loop
                else
                    allx = combnk(1:m.drugsTotal, j);
                end

                for i = 1:size(allx,1)
                    d = uint8(zeros(m.drugsTotal,1));
                    if(allx(i) ~= -1), d(allx(i,:)) = 1; end
                    v1 = m.calcObjVal(d,m.solutions.objs(1));
                    v2 = m.calcObjVal(d,m.solutions.objs(2));

                    withinconstraint = true;
                    %toxicity constraint
                    tox = m.dTox'*double(d);
                    for t = 1:length(tox)
                        if tox(t) > m.dToxMax(t), withinconstraint = false; end
                    end

                    if withinconstraint
                        m.solutions.all_cT = [m.solutions.all_cT; v1 v2];
                    else
                        m.solutions.all_cF = [m.solutions.all_cF; v1 v2];
                    end
                end
            end
        end
        
%% Single Objective Optimization
        function eff(m, idx)
            %MODIFIES: utopia, utopia_n, nadir, nadir_n
            cplex = Cplex(m.optModelName);
            cplex.Param.output.clonelog.Cur = 0;
            cplex.DisplayFunc = '';

            %objective
            eff = m.dEff*m.popstruct_e;
            lb = zeros(m.drugsTotal,1);
            ub = ones(m.drugsTotal,1);
            ctype = char(ones(1, m.drugsTotal)*('B'));
            cplex.addCols(eff, [], lb, ub, ctype);

            %fixed number of drugs
            if(m.drugsN > 0)
                s = ones(1,m.drugsTotal);
                cplex.addRows(m.drugsN, s, m.drugsN, 'drugsN');
            end

            %frontline drugs to include
            for c = 1:length(m.drugsInclude)
                s = zeros(1,m.drugsTotal);
                s(m.drugsInclude(c)) = 1;
                cplex.addRows(1, s, 1, 'drugsInclude');
            end

            %toxicity constraint
            for c = 1:size(m.dTox,2) %loop through each organ
                cplex.addRows(0, m.dTox(:,c)', m.dToxMax(c), 'DLT');
            end

            cplex.Model.sense = 'maximize';
            cplex.solve();
            m.solutions.utopia(idx) = cplex.Solution.objval;

            cplex.Model.sense = 'minimize';
            cplex.solve();
            m.solutions.nadir(idx) = cplex.Solution.objval;

            %normalize the points
            m.solutions.utopia_n(idx) = m.normalize(m.solutions.utopia(idx),idx);
            m.solutions.nadir_n(idx) = m.normalize(m.solutions.nadir(idx),idx);
        end
        
        function tox(m, idx)
            %MODIFIES: utopia, utopia_n, nadir, nadir_n
            cplex = Cplex(m.optModelName);
            cplex.Param.output.clonelog.Cur = 0;
            cplex.DisplayFunc = '';

            %objective
            tox = -m.dTox*ones(size(m.dTox,2),1);
            lb = zeros(m.drugsTotal,1);
            ub = ones(m.drugsTotal,1);
            ctype = char(ones(1, m.drugsTotal)*('B'));
            cplex.addCols(tox, [], lb, ub, ctype);

            %fixed number of drugs
            if(m.drugsN > 0)
                s = ones(1,m.drugsTotal);
                cplex.addRows(m.drugsN, s, m.drugsN, 'drugsN');
            end

            %frontline drugs to include
            for c = 1:length(m.drugsInclude)
                s = zeros(1,m.drugsTotal);
                s(m.drugsInclude(c)) = 1;
                cplex.addRows(1, s, 1, 'drugsInclude');
            end

            %toxicity constraint
            for c = 1:size(m.dTox,2) %loop through each organ
                cplex.addRows(0, m.dTox(:,c)', m.dToxMax(c), 'DLT');
            end

            cplex.Model.sense = 'maximize';
            cplex.solve();
            m.solutions.utopia(idx) = cplex.Solution.objval;

            cplex.Model.sense = 'minimize';
            cplex.solve();
            m.solutions.nadir(idx) = cplex.Solution.objval;

            %normalize the points
            m.solutions.utopia_n(idx) = m.normalize(m.solutions.utopia(idx),idx,true);
            m.solutions.nadir_n(idx) = m.normalize(m.solutions.nadir(idx),idx,true);
        end
        
%% Multi-Objective Optimization
        function d = efftox(m, rho, weights)
            %efficacy + frontline + toxicity
            %Augmented weighted Tchebycheff optimization
            %distance metric: Tchebycheff
            %INPUT: weights, size: 1 by 2
            %OUTPUT: drugs binary values, size: m.drugsTotal by 1; var type: uint8
            
            normfactor = m.solutions.normfactor;
            
            cplex = Cplex(m.optModelName);
            cplex.Param.output.clonelog.Cur = 0;
            cplex.Model.sense = 'minimize';
            cplex.DisplayFunc = '';

            %obj structure: [1:m.drugsTotal lambda]
            obj = zeros(m.drugsTotal+1,1);

            %lambda
            obj(end,1) = 1;

            %rho*L1 norm
            obj(1:end-1) = -rho*(m.dEff*m.popstruct_e - m.dTox*ones(size(m.dTox,2),1));

            lb = [zeros(m.drugsTotal,1); -inf];
            ub = [ones(m.drugsTotal,1); inf];
            ctype = char( [ones(1,m.drugsTotal)*('B') 'C'] );
            cplex.addCols(obj, [], lb, ub, ctype);

            %obj1 constraint
            eff = m.dEff*m.popstruct_e;
            eff = weights(1)*normfactor(1)*eff;
            s = [eff' 1];
            cplex.addRows(weights(1)*normfactor(1)*m.solutions.utopia(1), s, inf, 'efficacy_lambda');

            %obj2 constraint
            tox = -m.dTox*ones(size(m.dTox,2),1);
            tox = weights(2)*normfactor(2)*tox;
            s = [tox' 1];
            cplex.addRows(weights(2)*normfactor(2)*m.solutions.utopia(2), s, inf, 'toxicity_lambda');

            %fixed number of drugs
            if(m.drugsN > 0)
                s = ones(1,m.drugsTotal);
                cplex.addRows(m.drugsN, [s 0], m.drugsN, 'drugsN');
            end

            %frontline drugs to include
            for c = 1:length(m.drugsInclude)
                s = zeros(1,m.drugsTotal);
                s(m.drugsInclude(c)) = 1;
                cplex.addRows(1, [s 0], 1, 'drugsInclude');
            end

            %toxicity constraint
            for c = 1:size(m.dTox,2) %loop through each side effect
                cplex.addRows(-inf, [m.dTox(:,c)' 0], m.dToxMax(c), 'DLT');
            end

            cplex.solve();
            %cplex.writeModel([m.outputsdir m.optModelName '.lp']);

            %solutions are binary integers, save memory here by converting to uint8 (1 byte)
            d = uint8(cplex.Solution.x(1:m.drugsTotal));
        end
        
    end
end


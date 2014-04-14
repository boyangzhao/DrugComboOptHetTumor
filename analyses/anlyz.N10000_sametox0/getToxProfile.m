function toxprofile = getToxProfile(dTox, results)
% Tally the tox profile for the drug
% returns toxprofile, an array of struct; index for toxprofile corresponds
% to the number of drugs; each struct contains:
%    1) array of tox profile: each row corresponds to one solution;
%       column is arranged in same order as dTox (i.e. column 1 = Myelosuppression)
%    2) array of normal distri fit parameter values (empty), will be populated later when
%       generating hist in the next step. Row corresponds to toxicity.
%       Columns are [mu sigma].

    toxprofile = [];
    for i = 2:length(results) %will exclude result = 1, which presumably only has drugsN=0 solutions
        for j = 1:size(results(i).sols,1)

            tp = (dTox'*results(i).sols(j,:)')';

            drugsN = sum(results(i).sols(j,:),2);
            if drugsN == 0, continue; end
            
            while(drugsN > size(toxprofile,1))
                toxprofile = [toxprofile;struct('tp',[],'normfit',[])];
            end

            toxprofile(drugsN).tp = [toxprofile(drugsN).tp;tp];

        end
        if(mod(i,100.0) == 0), disp(num2str(i)); end
    end
end

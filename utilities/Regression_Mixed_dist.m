function [Statres, rEffects, Params] = Regression_Mixed_dist(PARRES,Time,Freq, indTotal)

nb  =1:numel(PARRES{1,1}.Model_reord);% number of boots to consider
VarNames = arrayfun(@(x) ['b' num2str(x)],nb,'uni',false);
ConNames = arrayfun(@(x) ['C' num2str(x)],1:numel(indTotal),'uni',false);% 30 connections

dists   = reshape(squeeze(mean(PARRES{1}.DistanceBoots.RFDists(:,:,:,:),3)),36,numel(PARRES{1,1}.Model_reord));
dists   = dists(indTotal,nb);
times = find(Time>-.2 );
%------------------------- reconstruct the data---------------------------
for cond = 1:numel(PARRES)
    
    Model_reord =   PARRES{cond}.Model_reord;
    for t = times
        t
        tic
        for f = 1:numel(Freq)
            for m = 1:numel(Model_reord)
                Data(:,:,m) = Model_reord{m}{1}.*Model_reord{m}{4}(f,:).*Model_reord{m}{5}(t,:);
                %Data(:,:,m) = Model_reord{m}{1}.*Model_reord{m}{5}(t,:);
            end
            %----------------------fit regression model--------------------
            for C = 1:4
                %(1) Format the data into table
                Data_table = array2table((squeeze(Data(:,C,nb))),'VariableNames',VarNames);
                %Data_table = array2table((squeeze(Data_con{cond}(:,C,nb))),'VariableNames',VarNames);
                Data_table = [table(ConNames','VariableNames',{'Conn'}) Data_table];
                Data_table.Conn = categorical(Data_table.Conn);
                Table_S = stack(Data_table,VarNames,...
                          'NewDataVariableName','iPDC',...
                          'IndexVariableName','Boots');

                % Distance data 
                Dist_table = array2table((squeeze(dists)),'VariableNames',VarNames);
                Dist_table = [table(ConNames','VariableNames',{'Conn'}) Dist_table];
                Dist_table.Conn = categorical(Dist_table.Conn);
                DTable_S = stack(Dist_table,VarNames,...
                          'NewDataVariableName','Dist',...
                          'IndexVariableName','Boots');

                TTable = join(Table_S,DTable_S,'keys',{'Conn','Boots'});

                TTable = TTable(~isnan(TTable.iPDC) & ~isnan(TTable.Dist),:);
                TTable = TTable(~isinf(TTable.iPDC) & ~isinf(TTable.Dist),:);

                %------------(2) Apply mixed-effect model with  random effect on intercept and slope
                lme_intercept_slope = fitlme(TTable,'iPDC ~ 1 + Dist + (1+Dist|Conn)');
                %lme_slope = fitlme(TTable,'iPDC ~ 1 + Dist +(Dist-1|Conn)');
                %Statres(cond,C,t,:,:)=[lme_intercept_slope.Coefficients.Estimate lme_intercept_slope.Coefficients.pValue];
                temp = lme_intercept_slope.Coefficients;
                temp2 = dataset2cell(temp);
                Statres{cond,C,t,f}= cell2mat(temp2(2:end,2:end));
                Params.Statres = temp2;
                
                % RANDOM EFFECTS
                [~,~,temp] = lme_intercept_slope.covarianceParameters;
                temp = temp{1};
                %[~,~,temp] = randomEffects(lme_intercept_slope);
                temp2 = dataset2cell(temp);
                rEffects{cond,C,t,f} =cell2mat(temp2(2:end,5:end));
                Params.rEffects = temp2;
                %-------------(3) Apply fixed-effect model---------------------------
                lm = fitlme(TTable,'iPDC ~ 1 + Dist');
                temp = lm.Coefficients;
                temp2 = dataset2cell(temp);
                Statlm{cond,C,t,f}=cell2mat(temp2(2:end,2:end));
                Params.Statlm = temp2;
                %--------------------(4) Compare models-----------------------------
                %Comp(cond,C,f,t,:) =  compare(lm, lme_intercept_slope).pValue;

            end
        end
        toc
    end
end


end
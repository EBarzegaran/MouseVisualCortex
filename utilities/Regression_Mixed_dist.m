function [Statres, rEffects, Params, Tables, RSquare] = Regression_Mixed_dist(PARRES,Time,Freq, indTotal,Formula, ConSel, KeepPercent, DODAI)

% 


if ~exist('KeepPercent','var')
    KeepPercent = 1;
end

if ~exist('DODAI','var')
    DODAI = true;
end

if DODAI==true
    KeepPercent = 1;
    ConSel = 1:numel(indTotal);
else
    if ~exist('ConSel','var')
        ConSel = 1:numel(indTotal);
    else
        indTotal = indTotal(ConSel);
    end
end

nb  =1:numel(PARRES{1,1}.Model_reord);% number of boots to consider
VarNames = arrayfun(@(x) ['b' num2str(x)],nb,'uni',false);
ConNames = arrayfun(@(x) ['C' num2str(x)],1:numel(indTotal),'uni',false);% 30 connections

dists   = reshape(squeeze(mean(PARRES{1}.DistanceBoots.RFDists(:,:,:,:),3)),36,numel(PARRES{1,1}.Model_reord))*9;% grid size 20
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
                if numel(Time)>1
                    TM = Model_reord{m}{5}(t,:).* mean(Model_reord{m}{2}).*mean(Model_reord{m}{3});
                else
                    TM = mean(Model_reord{m}{5}(:,:)).* mean(Model_reord{m}{2}).*mean(Model_reord{m}{3});
                end
                
                if numel(Freq)==1
                    Data(:,:,m) = Model_reord{m}{1}.*mean(Model_reord{m}{4}(f,:)).*TM;
                else
                    Data(:,:,m) = Model_reord{m}{1}.*Model_reord{m}{4}(f,:).*TM;
                end
                
            end
            %----------------------fit regression model--------------------
            for C = 1:4
                % remove 30% of connections with lowest connectivity values
                Data_temp = (squeeze(Data(ConSel,C,nb)));
                
                if ~DODAI
                    [~,Ind] = sort(sqrt(sum(Data_temp.^2,2)),'descend');
                    kpInd = Ind(1:round(numel(Ind)*KeepPercent));
                    Data_temp2 = (squeeze(Data(ConSel(kpInd),C,nb)));
                else
                    % Maybe later covert iPDC to DAI
                    temp = reshape(1:36,6,6);
                    a1 = tril(temp,-1);
                    a2 = tril(temp',-1);
                    DAI_I = [a1(:) a2(:)];
                    DAI_I(DAI_I(:,1)==0,:)=[];
                    [~,Loc1]=ismember(DAI_I(:,1),indTotal);
                    [~,Loc2]=ismember(DAI_I(:,2),indTotal);
                    DAI_I   = [Loc1 Loc2];
                    Data_temp2 = (Data_temp(DAI_I(:,1),:)-Data_temp(DAI_I(:,2),:));%./((Data_temp(DAI_I(:,1),:)+Data_temp(DAI_I(:,2),:)));
                    kpInd   = 1:numel(Loc1);
                end
                %(1) Format the data into table
                Data_table = array2table(Data_temp2,'VariableNames',VarNames);
                %Data_table = array2table((squeeze(Data_con{cond}(:,C,nb))),'VariableNames',VarNames);
                Data_table = [table(ConNames(kpInd)','VariableNames',{'Conn'}) Data_table];
                Data_table.Conn = categorical(Data_table.Conn);
                Table_S = stack(Data_table,VarNames,...
                          'NewDataVariableName','iPDC',...
                          'IndexVariableName','Boots');

                % Distance data 
                Dist_table = array2table((squeeze(dists(kpInd,:))),'VariableNames',VarNames);
                Dist_table = [table(ConNames(kpInd)','VariableNames',{'Conn'}) Dist_table];
                Dist_table.Conn = categorical(Dist_table.Conn);
                DTable_S = stack(Dist_table,VarNames,...
                          'NewDataVariableName','Dist',...
                          'IndexVariableName','Boots');

                TTable = join(Table_S,DTable_S,'keys',{'Conn','Boots'});

                TTable = TTable(~isnan(TTable.iPDC) & ~isnan(TTable.Dist),:);
                TTable = TTable(~isinf(TTable.iPDC) & ~isinf(TTable.Dist),:);
                
                Tables{cond,C}=TTable;

                %------------(2) Apply mixed-effect model with  random effect on intercept and slope
                lme_intercept_slope = fitlme(TTable,Formula);
                %lme_slope = fitlme(TTable,'iPDC ~ 1 + Dist +(Dist-1|Conn)');
                %Statres(cond,C,t,:,:)=[lme_intercept_slope.Coefficients.Estimate lme_intercept_slope.Coefficients.pValue];
                temp = lme_intercept_slope.Coefficients;
                temp2 = dataset2cell(temp);
                Statres{cond,C,t,f}= cell2mat(temp2(2:end,2:end));
                Params.Statres = temp2;
                
                % R-square
                RSquare{cond,C,t,f}=lme_intercept_slope.Rsquared;
                
                
                % RANDOM EFFECTS
%                 [A,B,temp] = lme_intercept_slope.covarianceParameters;
%                 temp = temp{1};
%                 temp2 = dataset2cell(temp);
                [~,~,temp2] = randomEffects(lme_intercept_slope);
                
                rEffects{cond,C,t,f} =temp2;%cell2mat(temp2(2:end,5:end));
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
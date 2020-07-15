function [R,P]= ConnectPDCCorr(PDC,RFDist,RFPval,P_th,diag,bsrm,time)


    NNode       = size(PDC,1);
    PDC = permute(PDC,[1 2 5 3 4]);
    PDC_temp    = reshape(PDC,[NNode*NNode*size(PDC,3) size(PDC,4) size(PDC,5)]);
    if bsrm
        PDC_temp = PDC_temp - mean(PDC_temp(:,:,(time<0) & (time>-.5)),3);
    end
    RFDist      = RFDist(:);
    
    RFPval = RFPval<P_th;
    for i = 1:size(RFPval,3)
        Temp = RFPval(1,:,i).*RFPval(2,:,i)';
        if diag
            Temp(1:length(Temp)+1:end)=0;
        end
        Select(:,:,i) = Temp==1;
    end
    Select = Select(:);
    
    for f = 1:size(PDC_temp,2)
        for t = 1:size(PDC_temp,3)
            [R(f,t) P(f,t)] = corr(PDC_temp(Select,f,t),RFDist(Select));
        end
    end
end
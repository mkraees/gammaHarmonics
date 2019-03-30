function [freqVals,psdST,psdBL,baseCorrectedLog10PSD] = getPSD(highRMSElectrodes,folderLFP,goodPos,stPos,blPos,mt)
% 3D Array - (Freq, trialIndex, electrodeInd)
    psdST = []; psdBL = []; baseCorrectedLog10PSD = [];
    for i = 1:length(highRMSElectrodes)
        load(fullfile(folderLFP,['elec' num2str(i)]),'analogData');
        stLFP = analogData(goodPos,stPos)';
        blLFP = analogData(goodPos,blPos)';
        [psdSTelec,freqVals] = mtspectrumc(stLFP,mt);
        psdBLelec = mtspectrumc(blLFP,mt);
        psdST = cat(3,psdST,psdSTelec);
        psdBL = cat(3,psdBL,psdBLelec);
        baseCorrectedLog10PSD = cat(3,baseCorrectedLog10PSD,log10(psdST)-log10(psdBL));
    end 
end
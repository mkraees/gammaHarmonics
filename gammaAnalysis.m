% function gammaAnalysis(subjectName, expType, stimType)

subjectName = 'tutu'; expType = 'Color'; stimType = '0';

    saveFolder = '/home/me/GammaHarmonicData/savedData';
    if isfile(fullfile(saveFolder,[subjectName expType stimType '.mat']))
        load(fullfile(saveFolder,[subjectName expType stimType '.mat']));
        disp('Variables loaded from Memory')
    else
        clearvars -except subjectName expType stimType saveFolder
        [highRMSElectrodes, goodPos, Fs, stPos, blPos, mt, folderLFP] = getSubjectDetails(subjectName, expType, stimType);
        [freqVals,psdST,psdBL,baseCorrectedLog10PSD, theta, RHO, stLFP, gammaSig, harmonicSig, gammaFreq, harmonicFreq, gammaAmp, harmonicAmp] = getPSDandPhase(highRMSElectrodes,folderLFP,goodPos,stPos,blPos,mt,Fs); %#ok<ASGLU>
        save(fullfile(saveFolder,[subjectName expType stimType '.mat']));
        disp('Variables Saved')
    end
    
    
% end


colormap jet;
pcolor(harmonicFreq); shading interp;
colorbar;
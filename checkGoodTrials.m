function [freqVals,psdST,psdBL,baseCorrectedLog10PSD, theta, RHO] = checkGoodTrials(subjectName, expType, stimType)

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
    
    I = divisors(length(goodPos));
    if length(I) <= 2
        I = divisors(length(goodPos)+1);
    end
    if mod(length(I),2)
        ROW = median(I);
        COL = ROW;
    else
        ROW = I(length(I)/2);
        COL = I(length(I)/2+1);
    end
    
    for i = 4
% %         f1 = figure('Name',['Electrode ' num2str(highRMSElectrodes(i))],'units','normalized','Position', [-1 0 1 1]);
        f2 = figure('Name',['Electrode ' num2str(highRMSElectrodes(i))],'units','normalized','Position', [1 0 1 1]);

        p = 1;
        for j = 1:length(goodPos)
%             
% %             figure(f1)
% %             subplot(ROW,COL,p)
% %             hold on
% %             plot(freqVals, baseCorrectedLog10PSD(:,j,i),'b','LineWidth',2);
% %             plot(freqVals, log10(psdST(:,j,i)),'r','LineWidth',1);
% %             plot(freqVals, log10(psdBL(:,j,i)),'k','LineWidth',1);
% %             xlim([0 150]);
% %             title(['Trial '  num2str(goodPos(j))]);
% %             hold off
%             
            [Mean, Std,~,~,~,Mode] = getPhaseProperties(RHO(:,j,i));
            Mean = Mean*180; Std = Std*180; Mode = Mode*180;
            
            figure(f2)
            subplot(ROW,COL,p)
            if (Mean >= 150 && Mean <= 210 && Std <=65 && Mode >= 150 && Mode <= 210)
                polarplot(theta',RHO(:,j,i),'r')
            else
                polarplot(theta',RHO(:,j,i))
            end
            
            title(['T-'  num2str(goodPos(j)) ' M-' num2str(round(Mean)) ' S-' num2str(round(Std)) ' Md-' num2str(round(Mode))]);
            p = p+1;
        end
        saveas(f2,[subjectName expType stimType 'PhaseDiff' '.tiff'])
    end
end


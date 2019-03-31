function figure2Gen(subjectName, expType, stimType)
% Generate Figure 2

    plotsFolder = '/home/me/GammaHarmonicData/Plots';
    makeDirectory(plotsFolder);
    saveFolder = '/home/me/GammaHarmonicData/savedData';
    makeDirectory(saveFolder);
    
    if isfile(fullfile(saveFolder,[subjectName expType stimType '.mat']))
        load(fullfile(saveFolder,[subjectName expType stimType '.mat'])); %#ok<LOAD>
        disp('Variables loaded from Memory')
    else
        clearvars -except subjectName expType stimType saveFolder plotsFolder
        [highRMSElectrodes, goodPos, Fs, stPos, blPos, mt, folderLFP, timeVals, stimPeriod, baselinePeriod, folderSourceString] = getSubjectDetails(subjectName, expType, stimType); %#ok<ASGLU>
        [freqVals,psdST,psdBL,baseCorrectedLog10PSD, theta, RHO, stLFP, gammaSig, harmonicSig, gammaFreq, harmonicFreq, gammaAmp, harmonicAmp] = getPSDandPhase(highRMSElectrodes,folderLFP,goodPos,stPos,blPos,mt,Fs); %#ok<ASGLU>
        save(fullfile(saveFolder,[subjectName expType stimType '.mat']));
        disp('Variables Saved')
    end

    G = find(gammaAmp >= median(median(gammaAmp)));
    H = find(harmonicAmp >= median(median(harmonicAmp)));

    GH = intersect(G,H);

    A = zeros(length(GH),2);
    for i = 1:length(GH)
        [trialInd, elecInd] = find(gammaAmp == gammaAmp(GH(i)));
        A(i,1) = elecInd;
        A(i,2) = trialInd;
    end

    selectedElecIndices = unique(A(:,1));
    selectedTrialIndices = unique(A(:,2));
    elecIndSelected = mode(A(:,1));
    trialIndSelected = mode(A(:,2));


    % [trialIndSelected, elecIndSelected] = find(harmonicAmp == max(max(harmonicAmp)));

    % PSD
    meanPSDBL = squeeze(mean(mean(psdBL,2),3));
    meanPSDST = squeeze(mean(mean(psdST,2),3));
    [freqGamma, freqHarmonic] = findGammaPeak(meanPSDST, freqVals);
    [singleTrialGamma, singleTrialHarmonic] = findGammaPeak(psdST(:,trialIndSelected,elecIndSelected), freqVals);

    figure('units','normalized','outerposition',[0 0 1 1])
    pos1 = [0.1 0.6 0.35 0.3];
    subplot('Position',pos1)

    plot(freqVals,log10(meanPSDBL)','k','LineWidth',1.5)
    hold on
    plot(freqVals,log10(meanPSDST)','r','LineWidth',1.5)
    plot(freqGamma,log10(meanPSDST(freqVals == freqGamma)),'ro','HandleVisibility','off');
    text(freqGamma,log10(meanPSDST(freqVals == freqGamma)),[' G = '  num2str(freqGamma) 'Hz']);
    plot(freqHarmonic,log10(meanPSDST(freqVals == freqHarmonic)),'ro','HandleVisibility','off');
    text(freqHarmonic,log10(meanPSDST(freqVals == freqHarmonic)),[' H = ' num2str(freqHarmonic) 'Hz']);

    plot(freqVals,log10(psdST(:,trialIndSelected,elecIndSelected))',':b','LineWidth',1.5)
    plot(singleTrialGamma,log10(psdST((find(freqVals == singleTrialGamma)),trialIndSelected,elecIndSelected)),'bo','HandleVisibility','off');
    plot(singleTrialHarmonic,log10(psdST((find(freqVals == singleTrialHarmonic)),trialIndSelected,elecIndSelected)),'bo','HandleVisibility','off')
    legend('BaseLine', 'Stimulus', 'SingleTrial', 'Location', 'Best')
    xlim([0 150])
    xlabel('Frequency (Hz)','FontSize', 12); ylabel('log10(PSD)','FontSize', 12);
    title([subjectName expType '-' stimType ' PSD'])
    set(gca, 'TickDir', 'out');
    hold off;


    % Raw LFP
    pos2 = [0.55 0.6 0.4 0.3];
    subplot('Position',pos2)
    stimTime = timeVals(stPos);
    hold on
    plot(stimTime, stLFP(:,trialIndSelected,elecIndSelected), 'b', 'LineWidth',1.5);
    plot(stimTime, gammaSig(:,trialIndSelected,elecIndSelected), 'g', 'LineWidth',1.5);
    plot(stimTime, harmonicSig(:,trialIndSelected,elecIndSelected), 'm', 'LineWidth',1.5);
    plot(stimTime, gammaSig(:,trialIndSelected,elecIndSelected) + harmonicSig(:,trialIndSelected,elecIndSelected), '--r', 'LineWidth',1.5);
    legend('LFP','G','H', 'G+H', 'Location', 'Best');
    xlabel('Time (s)','FontSize', 12)
    if strcmp(expType, 'Color')
        xlim([0.25 0.5])
    elseif  strcmp(expType, 'Length')
        xlim([0.5 0.75])
    end
    xlim([0.25 0.5])
    title('Single Trial LFP')
    set(gca, 'TickDir', 'out');
    hold off;

    % SingleElectrode SingleTrial
    subplot(234)
    rho = RHO(:,trialIndSelected,elecIndSelected);
    polarplot(theta',rho);
    title('Single Trial');

    % SingleElectrode AllTrials 
    subplot(235)
    polarplot(theta',mean(RHO(:,:,elecIndSelected),2));
    hold on
    polarplot(theta',mean(RHO(:,selectedTrialIndices,elecIndSelected),2));
    hold off
    title('SingleElectrode');

    %AllElectrodes AllTrials
    meanRHOall = mean(squeeze((mean(RHO(:,:,:),2))),2);
    meanRHOselected = mean(squeeze((mean(RHO(:,:,selectedElecIndices),2))),2);
    subplot(236)
    polarplot(theta',meanRHOall);
    hold on
    polarplot(theta',meanRHOselected);
    hold off
    title('AllElectrodes');
    
    saveas(gcf,fullfile(plotsFolder,[subjectName expType stimType 'Fig2' '.fig']))
    saveas(gcf,fullfile(plotsFolder,[subjectName expType stimType 'Fig2' '.tiff']))
    
end

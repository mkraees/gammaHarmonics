subjectName = 'tutu'; expType = 'Color'; stimType = 'Red';
% subjectName = 'alpa'; expType = 'Color'; stimType = 'Red';
% subjectName = 'alpa'; expType = 'Length'; stimType = 'con100';
% subjectName = 'kesari'; expType = 'Length'; stimType = 'con100';

[highRMSElectrodes, goodPos, Fs, stPos, blPos, mt, folderLFP, timeVals] = getSubjectDetails(subjectName, expType, stimType);

folderPair = fullfile('/home/me/GammaHarmonicData/ElecTrialPair');
PairStruct = load(fullfile(folderPair,[subjectName expType stimType 'GoodElectrodeTrialPair']));
trialPair = PairStruct.ElectrodeTrialPair;

elec = zeros(1,length(trialPair)); trialInd = zeros(1,length(trialPair));
for i = 1:length(trialPair)
    elec(i) = trialPair{i}(1);
    trialInd(i) = trialPair{i}(2);
end
electrodeInd = mode(elec);
electrodeNum = highRMSElectrodes(electrodeInd);
TrialInd = mode(trialInd);
TrialNum = goodPos(1);

subjectID = [subjectName expType];

% SingleTrial
i = electrodeNum;
load(fullfile(folderLFP,['elec' num2str(i)]),'analogData');
singleTrialLFP = analogData(TrialNum,stPos)';
[psdST, freqVals] = mtspectrumc(singleTrialLFP,mt); 
[singleTrialGamma, singleTrialHarmonic] = findGammaPeak(psdST, freqVals);

[freqVals, meanPSDST,meanPSDBL] = meanPSD(highRMSElectrodes, goodPos, folderLFP, stPos, blPos, mt);
[freqGamma, freqHarmonic] = findGammaPeak(meanPSDST, freqVals);

% PSD
figure('units','normalized','outerposition',[0 0 1 1])
pos1 = [0.1 0.6 0.35 0.3];
subplot('Position',pos1)
plot(freqVals,log10(meanPSDBL)','k','LineWidth',1.5)
hold on
plot(freqVals,log10(meanPSDST)','r','LineWidth',1.5)
plot(freqGamma,log10(meanPSDST(find(freqVals == freqGamma))),'ro','HandleVisibility','off');
text(freqGamma,log10(meanPSDST(find(freqVals == freqGamma))),[' G = '  num2str(freqGamma) 'Hz']);
plot(freqHarmonic,log10(meanPSDST(find(freqVals == freqHarmonic))),'ro','HandleVisibility','off');
text(freqHarmonic,log10(meanPSDST(find(freqVals == freqHarmonic))),[' H = ' num2str(freqHarmonic) 'Hz']);

plot(freqVals,log10(psdST)',':b','LineWidth',1.5)
plot(singleTrialGamma,log10(psdST(find(freqVals == singleTrialGamma))),'bo','HandleVisibility','off');
plot(singleTrialHarmonic,log10(psdST(find(freqVals == singleTrialHarmonic))),'bo','HandleVisibility','off')

legend('BaseLine', 'Stimulus', 'SingleTrial', 'Location', 'Best')
xlim([0 150])
xlabel('frequency (Hz)'); ylabel('log10(PSD)');
title([subjectID '-' stimType ' PSD'])
hold off;

% Raw LFP
pos2 = [0.55 0.6 0.4 0.3];
subplot('Position',pos2)
stimTime = timeVals(stPos);
[~, ~, gammaSig, harmonicSig] = getGammaPhaseDiff(singleTrialLFP, singleTrialGamma, singleTrialHarmonic, Fs);
hold on
plot(stimTime, singleTrialLFP, 'b', 'LineWidth',1.5);
plot(stimTime, gammaSig, 'g', 'LineWidth',1.5);
plot(stimTime, harmonicSig, 'm', 'LineWidth',1.5);
plot(stimTime, gammaSig + harmonicSig, '--r', 'LineWidth',1.5);
legend('SingleTrial','Gamma','Harmonic', 'Gamma + Harmonic', 'Location', 'Best');
xlabel('Time (s)')
if strcmp(expType, 'Color')
    xlim([0.25 0.5])
elseif  strcmp(expType, 'Length')
    xlim([0.5 0.75])
end
xlim([0.25 0.5])
title('Single Trial LFP')
hold off;

%%%%%%%%%%%%%%%%%%%%%%%

% SingleElectrode SingleTrial
trialIndex = find(goodPos == TrialNum);
[theta,Rho] = singleElectrodePhase(electrodeNum, goodPos, folderLFP, stPos, blPos, mt, Fs);
subplot(234)
rho = Rho(:,trialIndex)';
polarplot(theta,rho);
% legend('Trial');
title(['Elec = ' num2str(electrodeNum) ' Trial = ' num2str(TrialNum)]);


% SingleElectrode AllTrials 
selectSingleElec = 1;
[theta,Rho,GHpos] = singleElectrodePhase(electrodeNum, goodPos, folderLFP, stPos, blPos, mt, Fs, selectSingleElec);
subplot(235)
polarplot(theta',mean(Rho,2));
hold on
polarplot(theta',mean(Rho(:,GHpos),2));
% legend('All', 'Sel');
hold off
title(['Elec = ' num2str(electrodeNum)]);

        
%AllElectrodes AllTrials
RHO = [];
for i = 1:length(highRMSElectrodes)
    electrodeNum = i;
    dispFlag = 0;
    [theta,Rho] = singleElectrodePhase(electrodeNum, goodPos, folderLFP, stPos, blPos, mt, Fs);
    RHO = cat(2,RHO,mean(Rho,2));
end
meanRHO = mean(RHO,2); % mean of AllElec AllTrial
meanRHOselec = mean(RHO(:,elec),2); % mean of AllElec AllTrial
subplot(236)
polarplot(theta,meanRHO);
hold on
polarplot(theta,meanRHOselec);
% legend('All', 'Sel');
hold off
title(['AllElectrodes AllTrials']);
saveas(gcf,[subjectID stimType 'PhaseDiff' '.fig'])


% meanPSD
function [freqVals, meanPSDST,meanPSDBL] = meanPSD(highRMSElectrodes, goodPos, folderLFP, stPos, blPos, mt)
    meanPsdST = []; meanPsdBL = [];
    for i = 1:length(highRMSElectrodes)
        load(fullfile(folderLFP,['elec' num2str(i)]),'analogData');
        stLFP = analogData(goodPos,stPos)';
        blLFP = analogData(goodPos,blPos)';
        [psdST,freqVals] = mtspectrumc(stLFP,mt);
        psdBL = mtspectrumc(blLFP,mt);
        meanPsdST = cat(2, meanPsdST,mean(psdST,2));
        meanPsdBL = cat(2, meanPsdBL,mean(psdBL,2));
    end
    meanPSDST = mean(meanPsdST,2);
    meanPSDBL = mean(meanPsdBL,2);
end


% singleElectrodePhase
function [theta,Rho,GHpos] = singleElectrodePhase(electrodeNum, goodPos, folderLFP, stPos, blPos, mt, Fs, selectSingleElec)
    if ~exist('selectSingleElec', 'var'); selectSingleElec = 0; end
    i = electrodeNum;
    load(fullfile(folderLFP,['elec' num2str(i)]),'analogData');
    stLFP = analogData(goodPos,stPos)';
    blLFP = analogData(goodPos,blPos)';
    [psdST,freqVals] = mtspectrumc(stLFP,mt);
    psdBL = mtspectrumc(blLFP,mt);
    BLcorrectedPSD = log10(psdST)-log10(psdBL);
    
    Rho = []; peakGammaFreq = zeros(1,length(goodPos)); peakHarmonicFreq = zeros(1,length(goodPos));
    for j = 1:length(goodPos)
        [peakGammaFreq(i,j), peakHarmonicFreq(i,j)] = findGammaPeak(BLcorrectedPSD(:,j), freqVals);
        [theta, rho] = getGammaPhaseDiff(stLFP(:,j), peakGammaFreq(i,j), peakHarmonicFreq(i,j), Fs);
        Rho = cat(2,Rho,rho');
    end
    
    if selectSingleElec
         for j = 1:length(goodPos)
            [peakGammaFreq(j), peakHarmonicFreq(j),gammaAmp(j),harmonicAmp(j)] = findGammaPeak(BLcorrectedPSD(:,j), freqVals);
            [theta, rho] = getGammaPhaseDiff(stLFP(:,j), peakGammaFreq(j), peakHarmonicFreq(j), Fs);

            thetaPos = find(rho == max(rho));
            if mean(theta(thetaPos)) >=2.9 && mean(theta(thetaPos)) <= 3.4
                goodPhase(j) = 1;
            end
            Rho = cat(2,Rho,rho');
         end
        highGammaAmp = mean(gammaAmp);
        highHarmonicAamp = mean(harmonicAmp);

        highGammaPos = find(gammaAmp > highGammaAmp);
        highHarmonicPos = find(harmonicAmp > highHarmonicAamp);
        goodPhasePos = find(goodPhase == 1);

        GHpos0 = intersect(highHarmonicPos,highGammaPos);
        GHpos = intersect(GHpos0,goodPhasePos);
    end
end

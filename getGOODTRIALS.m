% function [Mean, Std] = getGoodElectodeTrialPair(subjectName, expType, stimType)
% subjectName = 'tutu'; expType = 'Color'; stimType = 'Red';
% subjectName = 'alpa'; expType = 'Color'; stimType = 'Red';
% subjectName = 'alpa'; expType = 'Length'; stimType = 'con100';
% subjectName = 'kesari'; expType = 'Length'; stimType = 'con100';


if ~exist('TW','var'); TW = 2; end
gridType = 'Microelectrode';

if strcmp(stimType, 'Red')
    GP = 1;
elseif strcmp(stimType, 'Green')
    GP = 120/10+1;
elseif strcmp(stimType, 'Blue')
    GP = 240/10+1;
elseif strcmp(stimType, 'con25')
    GP = 2;
elseif strcmp(stimType, 'con50')
    GP = 3;
elseif strcmp(stimType, 'con100')
    GP = 5;
end

subjectID = [subjectName expType];
if strcmp(subjectID,'alpaColor')
    subjectName = 'alpa';expDate = '301215'; protocolName = 'GRF_001'; % 488: Hue fullscreen
elseif strcmp(subjectID,'alpaLength')
    subjectName = 'alpa'; expDate = '120316'; protocolName = 'GRF_001';
elseif strcmp(subjectID,'kesariLength')
    subjectName = 'kesari'; expDate = '230716'; protocolName = 'GRF_002';
elseif strcmp(subjectID,'tutuColor')
    subjectName = 'tutu'; expDate = '191016'; protocolName = 'GRF_001'; % 111: Hue fullscreen
end

folderSourceString = fullfile('/home/me/GammaHarmonicData',expType);
folderBase = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
folderLFP = fullfile(folderBase,'segmentedData','LFP');
highRMSElectrodes = getHighRMSElectrodes(expType, subjectName, folderSourceString);
goodPosAll = getGoodPos(expType, folderBase);
goodPos = goodPosAll{GP};


if strcmp(expType, 'Color')
    stimPeriod = [0.25 0.75];% 500ms
    baselinePeriod = [-0.5 0]; %500 ms
elseif  strcmp(expType, 'Length')
    stimPeriod = [0.5 1.5]; %[0.5 2]
    baselinePeriod = [-1 0]; %[-1.5 0]
end

% TimeVals, FS and StimPos
load(fullfile(folderLFP,'lfpInfo'),'timeVals');
Fs = 1./(timeVals(2)-timeVals(1));
numPoints = round(diff(stimPeriod)*Fs);

stPos = find(timeVals>=stimPeriod(1),1) + (1:numPoints);
blPos = find(timeVals>=baselinePeriod(1),1) + (1:numPoints);

% Multi-taper parameters
tw = TW; fmax = 250;
mt.tapers = [tw (2*tw-1)];
mt.pad = -1; mt.Fs = Fs;
mt.fpass = [0 fmax];

gammaRangeHz = [30 70];
deltaF = 10; % Freqrange 20
n = 4; % order of butterWorth Filt
delta = deltaF;

[peakGammaFreq, peakHarmonicFreq, gammaAmp, harmonicAmp, goodPhase] = deal(zeros(length(highRMSElectrodes),length(goodPos)));
RsEsT = []; RsEmT = [];
for i = 1: length(highRMSElectrodes)
    load(fullfile(folderLFP,['elec' num2str(i)]),'analogData');
    stLFP = analogData(goodPos,stPos)';
    blLFP = analogData(goodPos,blPos)';
    [psdST,freqVals] = mtspectrumc(stLFP,mt);
    psdBL = mtspectrumc(blLFP,mt);
    BLcorrectedPSD = log(psdST)-log(psdBL);
    
    Rho = [];
    for j = 1:length(goodPos)
        [peakGammaFreq(i,j), peakHarmonicFreq(i,j),gammaAmp(i,j),harmonicAmp(i,j)] = findGammaPeak(BLcorrectedPSD(:,j), freqVals);
        [theta, rho] = GammaPhaseDiff(n, stLFP(:,j), peakGammaFreq(i,j), peakHarmonicFreq(i,j), delta, Fs);
        
        meanT = getPhaseProperties(rho);
        meanT = meanT*pi;
        if meanT >=2.9 && meanT <= 3.4
            goodPhase(i,j) = 1;
        end
        
        Rho = cat(2,Rho,rho');
    end
    RsEsT = cat(3,RsEsT,Rho); % SingleElec SingleTrial
    RsEmT = cat(2,RsEmT,mean(Rho,2)); % SingleElec mean of AllTrial
end

% Good Trials and Electrodes
highGammaAmp = max(mean(gammaAmp,2));
highHarmonicAamp = max(mean(harmonicAmp,2));

highGammaPos = find(gammaAmp > highGammaAmp);
highHarmonicPos = find(harmonicAmp > highHarmonicAamp);
goodPhasePos = find(goodPhase == 1);

GHpos0 = intersect(highHarmonicPos,highGammaPos);
% GHpos = GHpos0;
GHpos = intersect(GHpos0,goodPhasePos);

numElec = length(highRMSElectrodes);
ElectrodeTrialPair = cell(length(GHpos),1);
selectedRR = [];
for k = 1:length(GHpos)
    elec = mod(GHpos(k),numElec);
    trial = floor(GHpos(k)/numElec)+1;
    if elec == 0; elec = numElec; trial = trial-1; end
    ElectrodeTrialPair{k,1} = [elec trial];
    selectedRR = cat(2,selectedRR,mean(RsEsT(:,trial,elec),2));
end

save([subjectID stimType 'GoodElectrodeTrialPair'], 'ElectrodeTrialPair');

selectedRmEmT = mean(selectedRR,2); % mean of selected AllElec AllTrial

% Mean Phase plot
RmEmT = mean(RsEmT,2); % mean of AllElec AllTrial
Theta = theta';
[Mean, Std] = getPhaseProperties(RmEmT);

% figure()
% polarplot(Theta,RmEmT/sum(RmEmT))
% hold on
% polarplot(Theta,selectedRmEmT/sum(selectedRmEmT))
% title([subjectID '-' stimType ' - PhaseDiff - AllElec AllTrial' ' CircMean = ' num2str(Mean) ' CircStd = ' num2str(Std)])
% legend('All', 'Selected');
% hold off
% % saveas(gcf,[subjectID '-' stimType ' - PhaseDiff - AllElec AllTrial' '.fig'])

% end

% GammaPhaseDiff
function [theta, rho, gammaSig, harmonicSig] = GammaPhaseDiff(n, stLFP, peakGammaFreq, peakHarmonicFreq, delta, Fs)
    [B,A] = butter(n,[peakGammaFreq-delta, peakGammaFreq+delta]/(Fs/2));
    [D,C] = butter(n,[peakHarmonicFreq-delta, peakHarmonicFreq+delta]/(Fs/2));
    gammaSig = filtfilt(B,A,stLFP);
    harmonicSig = filtfilt(D,C,stLFP);

    G = hilbert(gammaSig);
    H = hilbert(harmonicSig);
    phaseDiff = (angle(H)-2*angle(G));
    [theta,rho] = rose(phaseDiff,25);
end

function [freqVals,psdST,psdBL,baseCorrectedLog10PSD, theta, RHO, stLFP, gammaSig, harmonicSig, gammaFreq, harmonicFreq, gammaAmp, harmonicAmp] = getPSDandPhase(highRMSElectrodes,folderLFP,goodPos,stPos,blPos,mt,Fs)
% PSD = 3D Array - (Freq, trialIndex, electrodeInd)
% RHO = 3D Array - (100, trialIndex, electrodeInd)
    n = 4; % order of butterWorth Filt
    delta = 10;
    
    gammaRangeHz =[30 70];
    
    psdST = []; psdBL = []; baseCorrectedLog10PSD = []; RHO = [];
    stLFP = []; gammaSig = []; harmonicSig = []; gammaFreq = []; harmonicFreq = []; gammaAmp = []; harmonicAmp = [];
    for i = 1:length(highRMSElectrodes)
        load(fullfile(folderLFP,['elec' num2str(i)]),'analogData');
        stLFPP = analogData(goodPos,stPos)';
        blLFP = analogData(goodPos,blPos)';
        [psdSTelec,freqVals] = mtspectrumc(stLFPP,mt);
        psdBLelec = mtspectrumc(blLFP,mt);
        psdST = cat(3,psdST,psdSTelec);
        psdBL = cat(3,psdBL,psdBLelec);
        baseCorrectedLog10PSD = cat(3,baseCorrectedLog10PSD,log10(psdST)-log10(psdBL));
        stLFP = cat(3,stLFP,stLFPP);
        
        gammaRangePos = intersect(find(freqVals>=gammaRangeHz(1)),find(freqVals<gammaRangeHz(2)));
        gammaRange = baseCorrectedLog10PSD(gammaRangePos,:,i);
        gammaAmpLog = max(gammaRange); % Amp difference from BL
        [gammaPosInd,~] = find(gammaRange==gammaAmpLog);
        gammaPos = gammaRangePos(gammaPosInd);
        peakGammaFreq = freqVals(gammaPos);
        peakHarmonicFreq = 2*peakGammaFreq; % Exact Harmonic
        gammaFreq = cat(2,gammaFreq,peakGammaFreq');
        harmonicFreq = cat(2,harmonicFreq,peakHarmonicFreq');
        harmonicPos = 2.*gammaPos-1;
        
        Rho = []; gammaSIG = []; harmonicSIG = []; gammaAmps = zeros(1,length(gammaPos)); harmonicAmps = zeros(1,length(gammaPos));
        for j =1:length(goodPos)
            [B,A] = butter(n,[peakGammaFreq(j)-delta, peakGammaFreq(j)+delta]/(Fs/2));
            [D,C] = butter(n,[peakHarmonicFreq(j)-delta, peakHarmonicFreq(j)+delta]/(Fs/2));
            
            gammaSignal = filtfilt(B,A,stLFPP(:,j));
            harmonicSignal = filtfilt(D,C,stLFPP(:,j));

            G = hilbert(gammaSignal);
            H = hilbert(harmonicSignal);
            phaseDiff = (angle(H)-2*angle(G));
            [theta,rho] = rose(phaseDiff,25); % 25*4
            Rho = cat(2,Rho,rho');
            
            gammaSIG = cat(2,gammaSIG,gammaSignal);
            harmonicSIG = cat(2,harmonicSIG,harmonicSignal);
            
            gammaAmps(j) = psdSTelec(gammaPos(j),j);
            harmonicAmps(j) = psdSTelec(harmonicPos(j),j);
        end
        RHO = cat(3,RHO,Rho);
        gammaSig = cat(3,gammaSig,gammaSIG);
        harmonicSig = cat(3,harmonicSig,harmonicSIG);
        gammaAmp = cat(2,gammaAmp,gammaAmps');
        harmonicAmp = cat(2,harmonicAmp,harmonicAmps');
    end 
end
function [peakGammaFreq, peakHarmonicFreq, gammaAmp, harmonicAmp, GHfreqRatio, GHampRatio, GHampRatioExact] = findGammaPeak(BLcorrectedPSD, freqVals, gammaRangeHz, deltaF)  
    if ~exist('gammaRangeHz','var'); gammaRangeHz =[30 70]; end
    if ~exist('deltaF','var'); deltaF = 10; end
    
    gammaRangePos = intersect(find(freqVals>=gammaRangeHz(1)),find(freqVals<gammaRangeHz(2)));
    gammaRange = BLcorrectedPSD(gammaRangePos);
    gammaAmp = max(gammaRange); % Amp difference from BL
    gammaPos = gammaRangePos(find(gammaRange==gammaAmp,1));
    peakGammaFreq = freqVals(gammaPos);
    
    harmonicEstFreq = 2*peakGammaFreq;
    peakHarmonicFreqExact = harmonicEstFreq;
    harmonicPosExact = freqVals == peakHarmonicFreqExact;
    harmonicAmpExact = BLcorrectedPSD(harmonicPosExact);
    
    harmonicRangePos = intersect(find(freqVals>=harmonicEstFreq-deltaF),find(freqVals<=harmonicEstFreq+deltaF));
    harmonicRange = BLcorrectedPSD(harmonicRangePos);
    harmonicAmp = max(harmonicRange); % Amp difference from BL
    harmonicPos = harmonicRangePos(find(harmonicRange==harmonicAmp,1));
    peakHarmonicFreq = freqVals(harmonicPos);
    
    GHfreqRatio = peakHarmonicFreq/peakGammaFreq; 
    GHampRatio = harmonicAmp/gammaAmp;
    GHampRatioExact = harmonicAmpExact/gammaAmp;
end
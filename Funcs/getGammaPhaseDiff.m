function [theta, rho, gammaSig, harmonicSig] = getGammaPhaseDiff(stLFP, peakGammaFreq, peakHarmonicFreq, Fs)
    n = 4; % order of butterWorth Filt
    delta = 10;
    [B,A] = butter(n,[peakGammaFreq-delta, peakGammaFreq+delta]/(Fs/2));
    [D,C] = butter(n,[peakHarmonicFreq-delta, peakHarmonicFreq+delta]/(Fs/2));
    gammaSig = filtfilt(B,A,stLFP);
    harmonicSig = filtfilt(D,C,stLFP);

    G = hilbert(gammaSig);
    H = hilbert(harmonicSig);
    phaseDiff = (angle(H)-2*angle(G));
    [theta,rho] = rose(phaseDiff,25);
end
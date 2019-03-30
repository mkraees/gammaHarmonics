function [highRMSElectrodes, goodPos, Fs, stPos, blPos, mt, folderLFP, timeVals, stimPeriod, baselinePeriod, folderSourceString] = getSubjectDetails(subjectName, expType, stimType) %#ok<STOUT>

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

    if isnan(str2double(stimType))
        if strcmp(stimType, 'Red')
            trialPos = 1;
        elseif strcmp(stimType, 'Green')
            trialPos = 120/10+1;
        elseif strcmp(stimType, 'Blue')
            trialPos = 240/10+1;
        elseif strcmp(stimType, 'con0')
            trialPos = 2;
        elseif strcmp(stimType, 'con25')
            trialPos = 2;
        elseif strcmp(stimType, 'con50')
            trialPos = 3;
        elseif strcmp(stimType, 'con75')
            trialPos = 4;
        elseif strcmp(stimType, 'con100')
            trialPos = 5;
        end
    else
        hue = str2double(stimType);
        trialPos = hue/10+1;
    end

    gridType = 'Microelectrode';
    folderSourceString = fullfile('/home/me/GammaHarmonicData',expType);
    folderBase = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
    folderLFP = fullfile(folderBase,'segmentedData','LFP');
    highRMSElectrodes = getHighRMSElectrodes(expType, subjectName, folderSourceString);
    goodPosAll = getGoodPos(expType, folderBase);
    goodPos = goodPosAll{trialPos};

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
    if ~exist('TW','var'); TW = 2; end
    tw = TW; fmax = 250;
    mt.tapers = [tw (2*tw-1)];
    mt.pad = -1; mt.Fs = Fs;
    mt.fpass = [0 fmax];

end
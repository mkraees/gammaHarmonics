function goodPos = getGoodPos(expType, folderBase)
    badTrialFile = fullfile(folderBase,'segmentedData','badTrials');
    load(badTrialFile,'badTrials');
    
    paramCombinationsFile = fullfile(folderBase,'extractedData','parameterCombinations');
    Params = load(paramCombinationsFile);
    cVals = Params.cValsUnique; oVals = Params.oValsUnique;

    if strcmp(expType, 'Color')
        goodPos = cell(1,length(oVals));
        c = find(cVals ==100); % contrast 100
        for o = 1:length(oVals)
            goodPos{o} = Params.parameterCombinations{1,1,1,1,o,c,1};
            goodPos{o} = setdiff(goodPos{o},badTrials);
        end
    elseif strcmp(expType, 'Length')
        goodPos = cell(1,length(cVals));
        for c = 1:length(cVals)
            goodPos{c} = Params.parameterCombinations{1,1,1,1,1,c,1};
            goodPos{c} = setdiff(goodPos{c},badTrials);
        end
    end
end
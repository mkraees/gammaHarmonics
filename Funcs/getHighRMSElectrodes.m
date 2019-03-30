function highRMSElectrodes = getHighRMSElectrodes(expType, subjectName, folderSourceString)
    if strcmp(expType, 'Color')
        folderElectrode = fullfile(folderSourceString,'analyzeElectrodes',subjectName);
        highRMSElectrodesStruct = load(fullfile(folderElectrode,'highRMSElectrodes'));
        highRMSElectrodes = highRMSElectrodesStruct.highRMSElectrodes;
    elseif  strcmp(expType, 'Length')
        folderRFData = fullfile(folderSourceString,'RFData');
        electrodeFile = fullfile(folderRFData, char([subjectName 'Microelectrode' 'RFData.mat']));
        electrodes = load(electrodeFile);
        highRMSElectrodes = electrodes.highRMSElectrodes;
    end
end
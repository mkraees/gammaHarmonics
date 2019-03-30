checkGoodTrials('tutu', 'Color', '0');
checkGoodTrials('alpa', 'Color', '0');

hue = 200:10:350;
Hue = num2str(hue(:));
tic
for i =1: length(hue)
    checkGoodTrials('tutu', 'Color', Hue(i,:));
    checkGoodTrials('alpa', 'Color', Hue(i,:));
end
toc

tic
checkGoodTrials('kesari', 'Length', 'con0');
checkGoodTrials('kesari', 'Length', 'con25');
checkGoodTrials('kesari', 'Length', 'con50');
checkGoodTrials('kesari', 'Length', 'con75');
checkGoodTrials('kesari', 'Length', 'con100');

checkGoodTrials('alpa', 'Length', 'con0');
checkGoodTrials('alpa', 'Length', 'con25');
checkGoodTrials('alpa', 'Length', 'con50');
checkGoodTrials('alpa', 'Length', 'con75');
checkGoodTrials('alpa', 'Length', 'con100');
toc
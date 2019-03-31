tic
figure2Gen('tutu', 'Color', '0');
figure2Gen('alpa', 'Color', '0');
hue = 10:10:90;
Hue = num2str(hue(:));
for i =1: length(hue)
    figure2Gen('tutu', 'Color', Hue(i,:));
    figure2Gen('alpa', 'Color', Hue(i,:));
    disp(hue(i));
end

close all
hue = 100:10:350;
Hue = num2str(hue(:));
for i =1: length(hue)
    figure2Gen('tutu', 'Color', Hue(i,:));
    figure2Gen('alpa', 'Color', Hue(i,:));
    disp(hue(i));
end
toc

close all
tic
disp('kesari');
figure2Gen('kesari', 'Length', 'con0');
figure2Gen('kesari', 'Length', 'con25');
figure2Gen('kesari', 'Length', 'con50');
figure2Gen('kesari', 'Length', 'con75');
figure2Gen('kesari', 'Length', 'con100');

disp('alpa')
figure2Gen('alpa', 'Length', 'con0');
figure2Gen('alpa', 'Length', 'con25');
figure2Gen('alpa', 'Length', 'con50');
figure2Gen('alpa', 'Length', 'con75');
figure2Gen('alpa', 'Length', 'con100');
toc

close all
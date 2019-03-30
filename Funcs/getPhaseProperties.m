function [Mean, Std, rayleighTest, C, Stat, Mode] = getPhaseProperties(R)
    C = []; 
    PP = round(R);
    for i = 1:100
        CC = zeros(1,PP(i));
        cc = 3.6*(i-1); % 360/100, 100 = 4*25;
        for j = 1:PP(i)
            CC(j) = cc;
        end
        C = cat(2,C,CC);
    end
    C = C/180;
    
    Stat = circ_stats(C);
    Mean = Stat.mean;
    Std = Stat.std;
    Mode = mode(C);
    rayleighTest = circ_rtest(C);  
end

sampletitle = cell(9,1);
% sampletitle{1} = 'Atto14';
sampletitle{1} = 'Atto15s7';
sampletitle{2} = 'T313ds7';
sampletitle{3} = 'T313ds8';
sampletitle{4} = 'T313ns7';
sampletitle{5} = 'T314ds7';
sampletitle{6} = 'T314ds8';
sampletitle{7} = 'ZF2W20s7';
sampletitle{8} = 'ZF2W21s7';
sampletitle{9} = 'ZF2W21s8';

elelist = cell(14,1);
elelist{1} = 'C';
elelist{2} = 'N';
elelist{3} = 'O';
elelist{4} = 'Na';
elelist{5} = 'Mg';
elelist{6} = 'P';
elelist{7} = 'S';
elelist{8} = 'Cl';
elelist{9} = 'K';
elelist{10} = 'Ca';
elelist{11} = 'Mn';
elelist{12} = 'Fe';
elelist{13} = 'Ni';
elelist{14} = 'Zn';

for q = 1:length(elelist)
    bigtotmass.(elelist{q}) = 0;
    bigmfrac.(elelist{q}) = 0;
end
bigtotmass.tot = 0;

for i = 1:9
    eval(['workingstruct = ' sampletitle{i} ';']);
    
    for k = 1:length(workingstruct.MixingOverview)
        for j = 1:length(elelist);
            
            bigtotmass.(elelist{j}) = workingstruct.MixingOverview(k).Mixing.totmass.(elelist{j}) + bigtotmass.(elelist{j});
        end
        bigtotmass.tot = workingstruct.MixingOverview(k).Mixing.totmass.tot + bigtotmass.tot;
        Dalist(k) = workingstruct.MixingOverview(k).Mixing.Da;
    end
    
    Hy = 0;
    for kk = 1:length(elelist)
        bigmfrac.(elelist{kk}) = bigtotmass.(elelist{kk}) ./ bigtotmass.tot;
        
        tempHy = - bigmfrac.(elelist{kk}) .* log(bigmfrac.(elelist{kk}));
        tempHy(isnan(tempHy) == 1) = 0;
        Hy = Hy + tempHy;
    end
    
    bigDa = mean(Dalist);
    bigDy = exp(Hy);
    
    bigChi = (bigDa - 1)./(bigDy - 1);
    
    workingstruct.bigDa = bigDa;
    workingstruct.bigDy = bigDy;
    workingstruct.bigChi = bigChi;
    
    eval([sampletitle{i} '= workingstruct;']);
    
end
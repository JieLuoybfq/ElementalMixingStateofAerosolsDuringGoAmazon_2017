function ulist = STXMEDXulist

ulist = cell(14,6);

%all mass absorption coefficients are in cm^2/g

%preC (278 eV)
ulist{1,1} = 1939; %C
ulist{2,1} = 3733; %N
ulist{3,1} = 5987; %O
ulist{4,1} = 17724; %Na
ulist{5,1} = 24856; %Mg
ulist{6,1} = 40933; %P
ulist{7,1} = 47416; %S
ulist{8,1} = 50255; %Cl
ulist{9,1} = 5699; %K
ulist{10,1} = 6904; %Ca
ulist{11,1} = 12943; %Mn
ulist{12,1} = 16235; %Fe
ulist{13,1} = 21200; %Ni
ulist{14,1} = 24509; %Zn

%postC (320 eV)
ulist{1,2} = 39697; %C
ulist{2,2} = 2545; %N
ulist{3,2} = 4210; %O
ulist{4,2} = 12513; %Na
ulist{5,2} = 17716; %Mg
ulist{6,2} = 30957; %P
ulist{7,2} = 36552; %S
ulist{8,2} = 39545; %Cl
ulist{9,2} = 47934; %K
ulist{10,2} = 5406; %Ca
ulist{11,2} = 9822; %Mn
ulist{12,2} = 12388; %Fe
ulist{13,2} = 16186; %Ni
ulist{14,2} = 19194; %Zn

%preN (398 eV)
ulist{1,3} = 24139; %C
ulist{2,3} = 1263; %N
ulist{3,3} = 2434; %O
ulist{4,3} = 7182; %Na
ulist{5,3} = 10384; %Mg
ulist{6,3} = 19922; %P
ulist{7,3} = 24165; %S
ulist{8,3} = 26711; %Cl
ulist{9,3} = 34403; %K
ulist{10,3} = 34887; %Ca
ulist{11,3} = 6327; %Mn
ulist{12,3} = 7643; %Fe
ulist{13,3} = 10384; %Ni
ulist{14,3} = 12692; %Zn

%postN (430 eV)
ulist{1,4} = 19930; %C
ulist{2,4} = 27505; %N
ulist{3,4} = 1997; %O
ulist{4,4} = 5878; %Na
ulist{5,4} = 8592; %Mg
ulist{6,4} = 16755; %P
ulist{7,4} = 20440; %S
ulist{8,4} = 22513; %Cl
ulist{9,4} = 29754; %K
ulist{10,4} = 32036; %Ca
ulist{11,4} = 5373; %Mn
ulist{12,4} = 6427; %Fe
ulist{13,4} = 8798; %Ni
ulist{14,4} = 10710; %Zn

%preO (525 eV)
ulist{1,5} = 12106; %C
ulist{2,5} = 17327; %N
ulist{3,5} = 1199; %O
ulist{4,5} = 3520; %Na
ulist{5,5} = 5116; %Mg
ulist{6,5} = 10528; %P
ulist{7,5} = 13002; %S
ulist{8,5} = 14112; %Cl
ulist{9,5} = 19346; %K
ulist{10,5} = 22012; %Ca
ulist{11,5} = 3469; %Mn
ulist{12,5} = 4092; %Fe
ulist{13,5} = 5636; %Ni
ulist{14,5} = 6648; %Zn

%postO (550 eV)
ulist{1,6} = 10763; %C
ulist{2,6} = 15558; %N
ulist{3,6} = 20708; %O
ulist{4,6} = 3126; %Na
ulist{5,6} = 4506; %Mg
ulist{6,6} = 9406; %P
ulist{7,6} = 11688; %S
ulist{8,6} = 12617; %Cl
ulist{9,6} = 17365; %K
ulist{10,6} = 19892; %Ca
ulist{11,6} = 3153; %Mn
ulist{12,6} = 3682; %Fe
ulist{13,6} = 5019; %Ni
ulist{14,6} = 5921; %Zn

end
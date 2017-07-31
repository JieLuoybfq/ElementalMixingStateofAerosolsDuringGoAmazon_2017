%[Mixing,Particles] = MixingState(Snew,filedir,names)
%
%Finding the mixing state of a set of STXM data that has already been run
%through the ParticleAnalysisGUI routine which includes LoadImageRawGui and
%CarbonMaps.
%
%Inputs
%------
%Snew (strucure output by CarbonMaps function)
%filedir string filedirectory found by uigetfile when ParticleAnalysis is is called and LoadFromImage is selected (maps)
%names filenames structure found by uigetfile when ParticleAnalysis is called and LoadFromImage is selected (maps)
%
%Outputs
%-------
%Mixing = Structure of mixing state parameters
%Mixing.TotalCompNum = total number of components present (only works w 3)
%Mixing.Hi = Shannon entropy for each particle given as a vector (1xN)
%Mixing.Ha = Average entropy per particle
%Mixing.Hy = Bulk entropy
%Mixing.Di = Per particle diversity
%Mixing.Da = Average diversity per particle
%Mixing.Dy = Average diversity of bulk
%Mixing.Db = interparticle diversity
%Mixing.MixStateChi = Entropy metric which quantifies the mixing state of a
%                       population.  Goes from 0 (internally mixed, all
%                       particles identically mixed) to 1 (externally mixed
%                       multiple particles, each a single component)
%
%Particles = Structure of particle information
%
%Non-Matlab native function dependencies
%------------
%LoadImageRawGui or LoadImageRawGuiMixingState
%AlignStack is compatible but not strictly dependent
%ODstack
%CarbonMaps
%Stacklabs is compatible but not strictly dependent
%ReadHdrWPixelSize
%
%Code by Matthew WF on 5/19/15
%
%
%Mixing state calculations taken from:
%Riemer, N., & West, M. (2013). Quantifying aerosol mixing state with entropy and diversity measures. Atmospheric Chemistry and Physics, 13(22), 11423-11439.

%maybe I should have a structure instead of a separate line for each calc?
%Thickness.soot etc or something?

function [Mixing,Particles] = MixingStateCNO(Snew)
% if abs(min(Snew.eVenergy) - 278) < 7 % 7 is an arbitrarily "close" value
%     if abs(max(Snew.eVenergy) - 320) < 7
%         errordlg('This is only the Carbon edge');
%     elseif abs(max(Snew.eVenergy) - 352.5) < 7
%         errordlg('This is the Carbon and Ca edge only');
%     elseif abs(max(Snew.eVenergy) - 430) < 7
%         errordlg('This includes the Carbon and Nitrogen edge only');
%     elseif abs(max(Snew.eVenergy) - 550) > 7
%         errordlg('This doesnt include the Oxygen postedge');
%     end
% else
%     errordlg('The Carbon pre-edge is missing');
% end



[~,idx278] = min(abs(Snew.eVenergy - 278));
[~,idx320] = min(abs(Snew.eVenergy - 320));
% % % [~,idx347] = min(abs(Snew.eVenergy - 347));
% % % [~,idx353] = min(abs(Snew.eVenergy - 352.5);
[~,idx400] = min(abs(Snew.eVenergy - 400));
[~,idx430] = min(abs(Snew.eVenergy - 430));
[~,idx525] = min(abs(Snew.eVenergy - 525));
[~,idx550] = min(abs(Snew.eVenergy - 550));

s278 = Snew.spectr(:,:,idx278); %spectra at listed energies
s320 = Snew.spectr(:,:,idx320);
s400 = Snew.spectr(:,:,idx400);
s430 = Snew.spectr(:,:,idx430);
s525 = Snew.spectr(:,:,idx525);
s550 = Snew.spectr(:,:,idx550);

u278 = 1939; %3554.7;
u320 = 39697; %17893;
u400 = 1243; %10739;
u430 = 27505; %12998;
u525 = 1199; %8022.1;
u550 = 20708; %14271;


D = 1.42; %density of alanine

position = Snew.position;
xstep = position.xstep;
ystep = position.ystep;
xstepcm = xstep ./ 10000; %converting to cm
ystepcm = ystep ./ 10000;
pixelarea = ystepcm .* xstepcm;
errpixelarea = sqrt((ystepcm.*ystepcm.*0.05).^2 + (xstepcm.*xstepcm.*0.05).^2);

cspectr = Snew.elemap.C; %s320-s278;
nspectr = Snew.elemap.N; %s430-s400;
ospectr = Snew.elemap.O; %s550-s525;
errcspectr = Snew.errmap.C;
errnspectr = Snew.errmap.N;
errospectr = Snew.errmap.O;

% cspectr(cspectr < 0) = 0;
% nspectr(nspectr < 0) = 0;
% ospectr(ospectr < 0) = 0;

cmu = u320-u278;
nmu = u430-u400;
omu = u550-u525;
errcmu = sqrt((u320.*0.05).^2 + (u278.*0.05).^2);
errnmu = sqrt((u430.*0.05).^2 + (u400.*0.05).^2);
erromu = sqrt((u550.*0.05).^2 + (u525.*0.05).^2);


mass.c = cspectr.*pixelarea./cmu;
mass.n = nspectr.*pixelarea./nmu;
mass.o = ospectr.*pixelarea./omu;
errmass.c = sqrt((pixelarea.*errcspectr./cmu).^2 + (cspectr.*errpixelarea./cmu).^2 + (-cspectr.*pixelarea.*errcmu./cmu.^2).^2);
errmass.n = sqrt((pixelarea.*errnspectr./nmu).^2 + (nspectr.*errpixelarea./nmu).^2 + (-nspectr.*pixelarea.*errnmu./nmu.^2).^2);
errmass.o = sqrt((pixelarea.*errospectr./omu).^2 + (ospectr.*errpixelarea./omu).^2 + (-ospectr.*pixelarea.*erromu./omu.^2).^2);


mass.tot = mass.c + mass.n + mass.o;
errmass.tot = sqrt((errmass.c.*mass.c).^2 + (errmass.n.*mass.n).^2 + (errmass.o.*mass.o).^2);

summass.c = sum(sum(mass.c));
summass.n = sum(sum(mass.n));
summass.o = sum(sum(mass.o));
summass.tot = sum(sum(mass.tot));

errsummass.c = sqrt(sum(sum(errmass.c.^2)));
errsummass.n = sqrt(sum(sum(errmass.n.^2)));
errsummass.o = sqrt(sum(sum(errmass.o.^2)));
errsummass.tot = sqrt(sum(sum(errmass.tot.^2)));

mfrac.c = summass.c ./ summass.tot;
mfrac.n = summass.n ./ summass.tot;
mfrac.o = summass.o ./ summass.tot;
errmfrac.c = sqrt((errsummass.c./summass.tot).^2 + (summass.c.*errsummass.tot./summass.tot.^2).^2);
errmfrac.n = sqrt((errsummass.n./summass.tot).^2 + (summass.n.*errsummass.tot./summass.tot.^2).^2);
errmfrac.o = sqrt((errsummass.o./summass.tot).^2 + (summass.o.*errsummass.tot./summass.tot.^2).^2);
Numparticles = max(max(Snew.LabelMat)); %total number of particles found

[row,column] = size(Snew.LabelMat); %row and column length of images, LabelMat chosen arbitrarily

Particles = struct('number', 1:Numparticles); %Defining structure
Particles.partmask = zeros(row,column,Numparticles); %preallocating matricies
Particles.partMmask.C = Particles.partmask;
Particles.partMmask.N = Particles.partmask;
Particles.partMmask.O = Particles.partmask;
Particles.errpartMmask.C = Particles.partmask;
Particles.errpartMmask.N = Particles.partmask;
Particles.errpartMmask.O = Particles.partmask;
Particles.Numparticles = Numparticles; %putting total number of particles into Particles struct
Particles.area = zeros(1,Numparticles);
Particles.errarea = Particles.area;
Particles.AED = Particles.area;
Particles.errAED = Particles.area;

%looping over number of particles and getting thickness of each component
%AND per particle
for i = 1:Numparticles
    partmask = Snew.LabelMat == i; %picking out each particle to mess with separately
    numpixels = sum(sum(partmask));
    Particles.area(1,i) = numpixels.*pixelarea;
    Particles.errarea(1,i) = numpixels.*errpixelarea;
    Particles.partmask(:,:,i) = partmask;
    Particles.partMmask.C(:,:,i) = mass.c .* partmask;
    Particles.partMmask.N(:,:,i) = mass.n .* partmask;
    Particles.partMmask.O(:,:,i) = mass.o .* partmask;
    Particles.errpartMmask.C(:,:,i) = errmass.c .* partmask;
    Particles.errpartMmask.N(:,:,i) = errmass.n .* partmask;
    Particles.errpartMmask.O(:,:,i) = errmass.o .* partmask;
    Particles.NumComp(i) = any(any(Particles.partMmask.C(:,:,i))) + any(any(Particles.partMmask.N(:,:,i))) + any(any(Particles.partMmask.O(:,:,i)));
end
Particles.AED = sqrt(4.*Particles.area./pi())*10000; %converting from cm to um  (10^6 / 10^2) = 10,000
Particles.errAED = sqrt((10000.*sqrt(4).*Particles.errarea./(2.*pi().*sqrt(Particles.area))).^2);

%masses of each component, a vector of masses indexed by particle number)
Particles.partM.C = permute(sum(sum(Particles.partMmask.C,1),2),[1,3,2]);
Particles.partM.N = permute(sum(sum(Particles.partMmask.N,1),2),[1,3,2]);
Particles.partM.O = permute(sum(sum(Particles.partMmask.O,1),2),[1,3,2]);
Particles.errpartM.C = permute(sum(sum(Particles.errpartMmask.C,1),2),[1,3,2]);
Particles.errpartM.N = permute(sum(sum(Particles.errpartMmask.N,1),2),[1,3,2]);
Particles.errpartM.O = permute(sum(sum(Particles.errpartMmask.O,1),2),[1,3,2]);

% % % 
% % % 
% % % 
% % % 
% Mixing State parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Particles.TotalM = Particles.partM.C + Particles.partM.N + Particles.partM.O;
Particles.errTotalM = sqrt((Particles.errpartM.C).^2 + (Particles.errpartM.N).^2 + (Particles.errpartM.O).^2);
FOVtotmass = sum(sum(mass.tot));
errFOVtotmass = sqrt(sum(sum(mass.tot.^2 .* errmass.tot.^2)));
Particles.Mfrac = Particles.TotalM ./ FOVtotmass;
Particles.errMfrac = sqrt((Particles.errTotalM./FOVtotmass).^2 + (-Particles.TotalM .* errFOVtotmass./FOVtotmass.^2).^2);

Particles.CompMfrac.C = Particles.partM.C ./ Particles.TotalM;
Particles.CompMfrac.N = Particles.partM.N ./ Particles.TotalM;
Particles.CompMfrac.O = Particles.partM.O ./ Particles.TotalM;
Particles.errCompMfrac.C = sqrt((Particles.errpartM.C./Particles.TotalM).^2 + (-Particles.partM.C .* Particles.errTotalM ./ Particles.TotalM.^2).^2);
Particles.errCompMfrac.N = sqrt((Particles.errpartM.N./Particles.TotalM).^2 + (-Particles.partM.N .* Particles.errTotalM ./ Particles.TotalM.^2).^2);
Particles.errCompMfrac.O = sqrt((Particles.errpartM.O./Particles.TotalM).^2 + (-Particles.partM.O .* Particles.errTotalM ./ Particles.TotalM.^2).^2);

%Defining total number of components
TotalCompNum = max(max(max(Particles.NumComp)));
Mixing = struct('TotalCompNum',TotalCompNum);

%Calculatin Shannon entropy per particle Hi
temp = zeros(3,Numparticles);
errtemp = zeros(3,Numparticles);
temp(1,:) = -Particles.CompMfrac.C .* log(Particles.CompMfrac.C);
temp(2,:) = -Particles.CompMfrac.N .* log(Particles.CompMfrac.N);
temp(3,:) = -Particles.CompMfrac.O .* log(Particles.CompMfrac.O);
errtemp(1,:) = sqrt(((-log(Particles.CompMfrac.C)-1).*Particles.errCompMfrac.C).^2);
errtemp(2,:) = sqrt(((-log(Particles.CompMfrac.N)-1).*Particles.errCompMfrac.N).^2);
errtemp(3,:) = sqrt(((-log(Particles.CompMfrac.O)-1).*Particles.errCompMfrac.O).^2);
temp(isnan(temp)) = 0;
Mixing.Hi = sum(temp);
Mixing.errHi = sqrt((errtemp(1,:)).^2 + (errtemp(2,:)).^2 + (errtemp(3,:)).^2);

%Calculating Shannon entropy for average particle Ha
weightedHi = Particles.Mfrac.*Mixing.Hi;
errweightedHi = sqrt((Particles.Mfrac.*Mixing.errHi).^2 + (Mixing.Hi.*Particles.errMfrac).^2);
Mixing.Ha = sum(weightedHi);
Mixing.errHa = sqrt(sum(errweightedHi.^2));

%Calculating Shannon entropy for bulk population Hy
Hy.C = -mfrac.c .* log(mfrac.c);
Hy.C(isnan(Hy.C)) = 0;
errHy.C = sqrt(((-log(mfrac.c)-1).*errmfrac.c).^2);

Hy.N = -mfrac.n .* log(mfrac.n);
Hy.N(isnan(Hy.N)) = 0;
errHy.N = sqrt(((-log(mfrac.n)-1).*errmfrac.n).^2);

Hy.O = -mfrac.o .* log(mfrac.o);
Hy.O(isnan(Hy.O)) = 0;
errHy.O = sqrt(((-log(mfrac.o)-1).*errmfrac.o).^2);
Mixing.Hy =  Hy.C + Hy.N + Hy.O;
Mixing.errHy = sqrt(errHy.C.^2 + errHy.N.^2 + errHy.O.^2);

%Calculating four diversity values: i - per particle diversity, a = average
%diversity, y = bulk population diversity, b = interparticle diversity
Mixing.Di = exp(Mixing.Hi);
Mixing.errDi = exp(Mixing.Hi).*Mixing.errHi;
Mixing.Da = exp(Mixing.Ha);
Mixing.errDa = exp(Mixing.Ha).*Mixing.errHa;
Mixing.Dy = exp(Mixing.Hy);
Mixing.errDy = exp(Mixing.Hy).*Mixing.errHy;
Mixing.Db = Mixing.Dy/Mixing.Da;
Mixing.errDb = sqrt((Mixing.errDy./Mixing.Da).^2 + (Mixing.Dy .* Mixing.errDa ./ Mixing.Da.^2).^2);

%Mixing State Index: 0 is externally mixed and 1 is internally mixed
Mixing.MixStateChi = (Mixing.Da - 1)/(Mixing.Dy - 1);
Mixing.errMixStateChi = sqrt((Mixing.errDa./(Mixing.Dy-1)).^2 + ((Mixing.Da-1).*Mixing.errDy./(Mixing.Dy-1).^2).^2);

%assignin('base','Mixing',Mixing)
%assignin('base','Particles',Particles)
end
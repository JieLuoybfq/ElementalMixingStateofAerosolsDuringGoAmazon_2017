%[totalmfrac] = mfractions(ParticlesOverview)
%
%Extracts mass fractions of soot, inorganics, and organics from the
%"ParticlesOverview structure into
%
%input
%-----
%ParticlesOverview
%ParticlesOverview.Particles.Numparticles
%ParticlesOverview.Particles.CompMfrac.C
%ParticlesOverview.Particles.CompMfrac.O
%ParticlesOverview.Particles.CompMfrac.N
%
%output
%------
%totalmfrac
%totalmfrac.soot
%totalmfrac.inorg
%totalmfrac.org
%
%non-matlab function synergies: MixingStateStats is meant to be run first

function [totalmfrac] = extractingmfracCNO(Dataset)

totalmfrac = struct('Carbon',0,'Nitrogen',0,'Oxygen',0);
totalnumberparticles = 0;
% names = fieldnames(Dataset.DataSet);
lstruc = length(Dataset);
for i = 1:lstruc
%     if isfield(Dataset.(names{i}),'Particles') == 0
%         continue
%     else
    totalnumberparticles = totalnumberparticles + Dataset(i).Particles.Numparticles;
    totalmfrac.Carbon = totalmfrac.Carbon + nansum(Dataset(i).Particles.CompMfrac.C);
    totalmfrac.errCarbon = sqrt(sum(Dataset(i).Particles.errCompMfrac.C.^2));
    totalmfrac.Nitrogen = totalmfrac.Nitrogen + nansum(Dataset(i).Particles.CompMfrac.N);
    totalmfrac.errNitrogen = sqrt(sum(Dataset(i).Particles.errCompMfrac.N.^2));
    totalmfrac.Oxygen = totalmfrac.Oxygen + nansum(Dataset(i).Particles.CompMfrac.O);
    totalmfrac.errOxygen = sqrt(sum(Dataset(i).Particles.errCompMfrac.O.^2));
%     end
end
totalmfrac.Carbon = totalmfrac.Carbon ./totalnumberparticles;
totalmfrac.errCarbon = totalmfrac.errCarbon ./totalnumberparticles;
totalmfrac.Nitrogen = totalmfrac.Nitrogen ./totalnumberparticles;
totalmfrac.errNitrogen = totalmfrac.errNitrogen ./totalnumberparticles;
totalmfrac.Oxygen = totalmfrac.Oxygen ./totalnumberparticles;
totalmfrac.errOxygen = totalmfrac.errOxygen ./totalnumberparticles;
end
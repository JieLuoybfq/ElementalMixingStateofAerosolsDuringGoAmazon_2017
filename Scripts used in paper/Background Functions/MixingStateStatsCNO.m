function [ Dataset ] = MixingStateStatsCNO(filedirs)
%SCRIPT
%Determination of simple statistics about mixing state and mass fractions
%Code by Matthew W. F. 5/19/15 @ University of the Pacific

%non-matlab function dependencies
%--------------------------------
%uipickfiles (internet)
%LoadImageRawGuiMixingStateOutput (file by RCM modified by MWF)
%mfractions (file by MWF)


ldirs = length(filedirs);
dirnames = cell(ldirs,1);
emptycell = cell(ldirs,1);
for i = 1:ldirs
    [~,dirnames{i},~] = fileparts(filedirs{i});
    dirnames{i} = ['FOV' dirnames{i}];
    idx = strfind(dirnames{i},'-'); %cell2struct doesn't like hyphens
    dirnames{i}(idx) = '_';
end

Dataset = cell2struct(emptycell,dirnames,1);
ParticlesOverview = zeros(1,ldirs);
MixingOverview = zeros(1,ldirs);


% filedirs = uipickfiles; %ui window for picking multiple directories
%PREALLOCATE MIXINGOVERVIEW?
MixingOverview = struct('DataSet',cell(1,length(filedirs)));
ParticlesOverview = struct('DataSet',cell(1,length(filedirs)));
tempDa = zeros(1,length(filedirs));
tempDy = zeros(1,length(filedirs));
tempDb = zeros(1,length(filedirs));
tempCHI = zeros(1,length(filedirs));

if any(exist('sillystring','file'))
	hwait = waitbar(0,sillystring);
else
	hwait = waitbar(0,'plz w8');
end

for i = 1:length(filedirs); %looping through each selected directory
    tempfiledir = strcat(filedirs{i},'\');
    cd(filedirs{i}); %moving to each directory
    tempfilenames = ls; %listing out file names
    cnt = 1;
    hdrcnt = 0;
    for j = 1:size(tempfilenames,1) %looping through each file name and picking out .hdr and .xim files ONLY 
                                    %(this allows for lots of other crap in the folder) and then building 
                                    %the FileNames cell array
        if any(strfind(tempfilenames(j,:),'.hdr'))==1
           FileNames{1,cnt} = strtrim(tempfilenames(j,:));    %i'm not sure how to preallocate here without another if loop, might not be faster
           cnt = cnt + 1;
           hdrcnt = hdrcnt + 1;
        elseif any(strfind(tempfilenames(j,:),'.xim'))==1
            FileNames{1,cnt} = strtrim(tempfilenames(j,:));
            cnt = cnt + 1;
        end
    end
    if hdrcnt > 1 %S and Snew hidden, maybe make an Snew overview struct?
        [S,Snew,Mixing,Particles] = LoadImageRawMixingStateOutputCNO(tempfiledir,FileNames);
    elseif hdrcnt == 1
        [S,Snew,Mixing,Particles] = SingStackProcMixingStateOutputCNO(tempfiledir);
    end
    
    Dataset.(dirnames{i}).S = S;
	Dataset.(dirnames{i}).Snew = Snew;
	Dataset.(dirnames{i}).Mixing = Mixing;
	Dataset.(dirnames{i}).Particles = Particles;
	Dataset.(dirnames{i}).Directory = filedirs{i};

    
    MixingOverview(1,i).DataSet = tempfiledir;
    MixingOverview(1,i).Mixing = Mixing;
    MixingOverview(1,i).Mixing.Numparticles = length(Mixing.Di);
    tempDa(i) = MixingOverview(1,i).Mixing.Da;
    tempDy(i) = MixingOverview(1,i).Mixing.Dy;
    tempDb(i) = MixingOverview(1,i).Mixing.Db;
    tempCHI(i) = MixingOverview(1,i).Mixing.MixStateChi;
    ParticlesOverview(1,i).DataSet = tempfiledir;
    ParticlesOverview(1,i).Particles = Particles;
    clear FileNames
    waitbar(i/length(filedirs));
end
close(hwait);


DaStats = SimpleStats(tempDa);
DyStats = SimpleStats(tempDy);
DbStats = SimpleStats(tempDb);
ChiStats = SimpleStats(tempCHI);

MixStateStats = struct('DaStats',DaStats,'DyStats',DyStats,'DbStats',DbStats,'ChiStats',ChiStats);
Dataset.MixStateStats = MixStateStats;


totalmfrac = extractingmfracCNO(ParticlesOverview);
Dataset.totalmfrac = totalmfrac;


% assignin('base','MixingOverview',MixingOverview);
% assignin('base','MixStateStats',MixStateStats);
% assignin('base','ParticlesOverview',ParticlesOverview);

% clear i j tempCHI tempDa tempDb tempDy hdrcnt cnt Mixing Particles filedirs;
% clear S Snew tempfiledir tempfilenames ChiStats DaStats DbStats DyStats;

end
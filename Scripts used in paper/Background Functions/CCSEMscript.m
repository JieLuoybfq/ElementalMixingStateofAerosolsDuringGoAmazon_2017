function [semjob] = CCSEMscript(Optionalfiledir)
%[semjob] = CCSEMscript()
%semjob is a struct, CCSEMscript is called without input arguments
%
%navigate to the folder containing all data (all stubs) related to a single
%CCSEM job and click "add"

%Outputs (semjob is replaced by whatever user chooses output variable name
%           to be)
%===============================================================
%semjob.stub#.compositepic = all fields of a stub pieced together
%semjob.stub#.fieldposinfo.xstgnum = number of columns
%semjob.stub#.fieldposinfo.ystgnum = number of rows
%semjob.stub#.fieldposinfo.fieldnum = total number of fields
%semjob.stub#.imgdata = a 3d matrix which is (m X n X fieldnum)
%                           where m X n is the img of 1 field
%semjob.stub#.semdata = a large per-particle matrix of edx data
%semjob.stub#.colheaders = the column headings of semdata
%semjob.stub#.fieldend = a list of the last particle index of each field (inclusive)
%semjob.stub#.numparticles = a list of the number of particles in each field

%Coded by Matthew W. Fraund 2/3/16 from the University of the Pacific in
%Stockton, CA.


%this checks for the existence of the online matlab script "uipickfiles"
%which is great and you should download if you don't have it.  It uses
%matlab's native uigetdir if you don't have it.

if nargin == 0
    
    cd('C:\Users');
    userdir = dir;
    userdir = userdir(3:end);
    for i = 1:length(userdir)
        if strcmp(userdir(i).name,'Katy-Ann') == 1
            cd('C:/Users/Katy-Ann/Google Drive/Projects/GoAmazon/SEM_EDX_Data/Matthew_EDX');
            break
        elseif strcmp(userdir(i).name,'Emily') == 1
            cd('C:/Users/Emily/Google Drive/Projects/GoAmazon/SEM_EDX_Data/Matthew_EDX');
            break
        elseif strcmp(userdir(i).name,'mulecenter78') == 1
            cd('C:/Users/mulecenter78/Google Drive/Projects/GoAmazon/SEM_EDX_Data/Matthew_EDX');
            break
        elseif strcmp(userdir(i).name,'Shodan') == 1
            cd('C:/Users/Shodan/Google Drive/Projects/GoAmazon/SEM_EDX_Data/Matthew_EDX');
            break
        end
    end
    
    if any(exist('uipickfiles','file'))
        ccdir = uipickfiles('Prompt',...
            'Pick directory containing all files from 1 SEM job');
        cd(ccdir{1});
    else
        ccdir = uigetdir;
        cd(ccdir);
    end
    
elseif nargin == 1
    ccdir = Optionalfiledir;
    cd(ccdir{1});
end
[~,foldername,~] = fileparts(ccdir{1});

%%%temporary for optimizing SURF matching
% % ccdir = 'C:\Users\Katy-Ann\Google Drive\SEMTEST\CCSEMof150913001+2';
% % cd(ccdir);
% % [~,foldername,~] = fileparts(ccdir);
%%%


%finding the csv file (data file)
jobls = ls;
stubnum = 0;
for i = 1:size(jobls,1)
    stubidx = strfind(strtrim(jobls(i,:)),'stub');
    if isdir(strtrim(jobls(i,:))) && ~isempty(stubidx)
        stubnum = stubnum + 1;
        dirrow(stubnum) = i;
    end
end

for kk = 1:stubnum %this loops over multiple stubs but I haven't really used it, so attempting multi-stub analysis is probably buggy =(
    cd(strtrim(jobls(dirrow(kk),:)));
    stubls = ls;
    dirnum = 1;
    fielddirlist = cell(1,1);
    for i = 1:size(stubls,1)
        idx = strfind(strtrim(stubls(i,:)),'.csv');
        idx2 = strfind(strtrim(stubls(i,:)),'.txt');
        flddiridx = strfind(strtrim(stubls(i,:)),'fld');
        if ~isempty(idx)
            csvrow = i;
        elseif ~isempty(idx2)
            txtrow = i;
        elseif ~isempty(flddiridx)
            fielddirlist{dirnum,1} = strtrim(stubls(i,:));
            dirnum = dirnum + 1;
        end
    end
    
    csvfile = strtrim(stubls(csvrow,:));
    txtfile = strtrim(stubls(txtrow,:));
    
% % %     %importing data file while removing text lines at top
% % %     DELIMITER = ',';
% % %     HEADERLINES = 15; %15 lines for stub csv, 13 lines for job csv file
% % %     
% % %     % Import the file
% % %     newData1 = importdata(csvfile, DELIMITER, HEADERLINES);
% % %     
% % %     % Create new variables in the base workspace from those fields.
% % %     vars = fieldnames(newData1);
% % %     for i = 1:length(vars) %this makes 3 variables: data, colheaders, textdata
% % %         eval(sprintf('%s',vars{i},' = newData1.(vars{i});'));
% % %     end

    [data,colheaders,textdata] = readCCSEMcsv(csvfile);
    [pixelsize] = readCCSEMstubsummary(txtfile);
    
    fieldnum = max(data(:,2)); %data(:,2) is the "field #" column
    
    %because fields are scanned from right to left before changing
    %y direction, this xstgnum is the total number of fields in this
    %scan and the intuitive xstgnum is calculated later.  
    %An exception exists for 1xn scans, where the x stage value never changes
    %which is dealt with a bit later
    xstages = data(1,4);
    xstgnum = 1;
    fldidx = 1;
    fieldend = [];
    for i = 1:size(data,1)
        if data(i,4) == xstages(xstgnum) %if it's the same stage, do nothing
            continue
        else
            fieldend(fldidx) = i-1;
            fldidx = fldidx + 1;
            xstages = [xstages data(i,4)];
            xstgnum = xstgnum + 1;
        end
    end
    fieldend(fldidx) = i; %adding the index for the final entry
    
    %y dimension only changes after all x translation is done for
    %that given y value.  ystgnum is how many rows the scan job contained
    ystages = data(1,5);
    ystgnum = 1;
    
    if xstgnum == 1 %special loop counting fldidx/fieldend if xstgnum never changed in loop above (i.e. 1xn stubs)
        for i = 1:size(data,1)
            if data(i,5) == ystages(ystgnum)
                continue
            else
                fieldend(fldidx) = i-1;
                fldidx = fldidx + 1;
                ystages = [ystages data(i,5)];
                ystgnum = ystgnum + 1;
            end
        end
        fieldend(fldidx) = i; %index for final entry
    else %if the number of x stages is >1
        for i = 1:size(data,1)
            if data(i,5) == ystages(ystgnum)
                continue
            else
                ystages = [ystages data(i,5)];
                ystgnum = ystgnum + 1;
            end
        end
    end
    
    if xstgnum == 1
    else
        xstgnum = fieldnum./ystgnum; %necessary to deal with how xstgnum is calculated above
    end
    
    
    %counting particles in a given field
    numparticles = zeros(fieldnum,1);
    lastfieldidx = 0;
    for i = 1:length(fieldend);
        numparticles(i,1) = fieldend(i) - lastfieldidx;
        lastfieldidx = fieldend(i);
    end
    
    imgdata = []; %needed for first cat command to work
    for k = 1:fieldnum %stitching together SEM pictures
        
        %Imports .tif files
        fileToRead1 = sprintf('%s',fielddirlist{k,1},'/','search.tif');
        
        % Import the file
        newData2 = importdata(fileToRead1);
        
        % Create new variables in the base workspace from those fields.
        vars = fieldnames(newData2);
        
        %this makes 2 variables: cdata and colormap.
        %    cdata is important, colormap isnt
        for i = 1:length(vars)
            eval(sprintf('%s',vars{i},'= newData2.(vars{i});'));
        end
        
        %each tif img, is stored in individual 2d matricies
        imgdata = cat(3,imgdata,cdata);
    end
    
    compositepic = []; %needed for cat command to work for the first matrix
    for j = 1:ystgnum
        rowdata = [];
        rownums = ((xstgnum*(j-1)+1):xstgnum*j); %this will create a vector of numbers corresponding to a single row (1 2 3) or (4 5 6) etc
        for q = 1:length(rownums)
            rowdata = cat(2,rowdata,imgdata(:,:,rownums(q)));
        end
        compositepic = cat(1,compositepic,rowdata);
    end
    
    errpixelsize = pixelsize.*0.05; %2 um repeatability
    
    xsize = size(compositepic,2) .* pixelsize;
    errxsize = size(compositepic,2).*errpixelsize;
    ysize = size(compositepic,1) .* pixelsize;
    errysize = size(compositepic,1).*errpixelsize;
    
    stubname = sprintf('%s','stub',num2str(stubnum)); %this is me trying to code for multiple stubs
    
    
    fieldposinfo.xstgnum = xstgnum;
    fieldposinfo.ystgnum = ystgnum;
    fieldposinfo.fieldnum = fieldnum;
    
    eval(sprintf('%s',stubname,'.fieldposinfo = fieldposinfo;'));
    eval(sprintf('%s',stubname,'.compositepic = compositepic;'));
    eval(sprintf('%s',stubname,'.imgdata = imgdata;'));
    eval(sprintf('%s',stubname,'.semdata = data;'));
    eval(sprintf('%s',stubname,'.colheaders = colheaders;'));
    eval(sprintf('%s',stubname,'.fieldend = fieldend;'));
    eval(sprintf('%s',stubname,'.numparticles = numparticles;'));
    eval(sprintf('%s',stubname,'.pixelsize = pixelsize;'));
    eval(sprintf('%s',stubname,'.xsize = xsize;'));
    eval(sprintf('%s',stubname,'.ysize = ysize;'));
    eval(sprintf('%s',stubname,'.errxsize = errxsize;'));
    eval(sprintf('%s',stubname,'.errysize = errysize;'));
    eval(sprintf('%s',stubname,'.errpixelsize = errpixelsize;'));
    eval(sprintf('%s',stubname,'.foldername = foldername;'));
    
    semjob.(stubname) = eval(sprintf('%s',stubname)); %outputing a single structure is much cleaner than multiple variables
end

end

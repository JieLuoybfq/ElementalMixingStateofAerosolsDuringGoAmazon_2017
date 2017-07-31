function matchingsemstxm()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      SEM Initialization   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%finding SEM directory
semdata = CCSEMscript;
sempixsize = semdata.stub1.pixelsize;
if length(fieldnames(semdata)) > 2;
    errordlg('Too many SEM jobs chosen, only pick 1');
end


%changing the pixel positions in entire data set into row/column positions
coldata = semdata.stub1.semdata(:,6);
rowdata = (size(semdata.stub1.imgdata,1) - semdata.stub1.semdata(:,7) + 1);

semdata.stub1.semdata(:,6) = rowdata; %NOTE in colheaders, column 6 is x_cent and 7 is y_cent
semdata.stub1.semdata(:,7) = coldata; %after this transformation, column 6 is rows and 7 is cols

semdata.stub1.colheaders{6} = 'row_cent';
semdata.stub1.colheaders{7} = 'col_cent';



%converting particle x and y PIXEL values with 0,0 based on field number into x and
%y PIXEL values with 0,0 based on the entire composite image
numxfields = semdata.stub1.fieldposinfo.xstgnum;
numyfields = semdata.stub1.fieldposinfo.ystgnum;
fieldxsize = size(semdata.stub1.imgdata,2); %size(...,2) is cols, which is x data
fieldysize = size(semdata.stub1.imgdata,1); %rows is y data
for ii = 1:size(semdata.stub1.semdata,1);
    curfield = semdata.stub1.semdata(ii,2); %this is field number independent of how fields are split horizontally/vertically
    yfields = 1;
    tempfield = curfield;
    while tempfield > numxfields
        tempfield = tempfield - numxfields;
        yfields = yfields + 1;
    end
    xfields = tempfield;
    
    semdata.stub1.semdata(ii,6) = semdata.stub1.semdata(ii,6) + fieldysize.*(yfields-1);
    semdata.stub1.semdata(ii,7) = semdata.stub1.semdata(ii,7) + fieldxsize.*(xfields-1);
end

%a little image enhancement
SEMimageraw = mat2gray(semdata.stub1.compositepic);
SEMimage = SEMimageraw;
SEMimage = imadjust(SEMimage,[0;1],[0;1],1);


% % %changing the pixel positions in entire data set into row/column positions
% % coldata = semdata.stub1.semdata(:,6);
% % rowdata = (size(SEMimage,1) - semdata.stub1.semdata(:,7) + 1);
% %
% % semdata.stub1.semdata(:,6) = rowdata; %NOTE in colheaders, column 6 is x_cent and 7 is y_cent
% % semdata.stub1.semdata(:,7) = coldata; %after this transformation, column 6 is rows and 7 is cols
% %
% % semdata.stub1.colheaders{6} = 'row_cent';
% % semdata.stub1.colheaders{7} = 'col_cent';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%        STXM Initilization     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finding STXM directory
% stxmdir = uipickfiles('Prompt','Select Folder Containing STXM Data');
% cd('C:/Users/Katy-Ann/Google Drive/Projects/GoAmazon/All Amazon STXM Data');

cd('C:\Users');
userdir = dir;
userdir = userdir(3:end);
for i = 1:length(userdir)
    if strcmp(userdir(i).name,'Katy-Ann') == 1
        cd('C:/Users/Katy-Ann/Google Drive/Projects/GoAmazon/All Amazon STXM Data');
        break
    elseif strcmp(userdir(i).name,'Emily') == 1
        cd('C:/Users/Emily/Google Drive/Projects/GoAmazon/All Amazon STXM Data');
        break
    elseif strcmp(userdir(i).name,'mulecenter78') == 1
        cd('C:/Users/mulecenter78/Google Drive/Projects/GoAmazon/All Amazon STXM Data');
        break
    elseif strcmp(userdir(i).name,'Shodan') == 1
        cd('C:/Users/Shodan/Google Drive/Projects/GoAmazon/All Amazon STXM Data');
        break
    end
end


if any(exist('uipickfiles','file'))
    stxmdir = uipickfiles('Prompt',...
        'Select Folder Containing STXM Data');
    cd(stxmdir{1});
else
    stxmdir = uigetdir;
    blank = cell(1); %this is to make sure it's in the same format as uipickfiles outputs
    blank{1} = stxmdir;
    stxmdir = blank;
    cd(stxmdir);
end

%%%%%%temporary for optimizing surf stuff
% % stxmdir = cell(1);
% % stxmdir{1} = 'C:\Users\Katy-Ann\Google Drive\SEMTEST\ShortStack532_150913002';
%%%%%%

if length(stxmdir) > 2;
    errordlg('Too many STXM FOV"s chosen, only pick 1');
end

[~,stxmfolder,~] = fileparts(stxmdir{1});
stxmdata = MixingStateStatsCNO(stxmdir);
stxmfieldnames = fieldnames(stxmdata);

%pixel size of stxm image
stxmpixsize = stxmdata.(stxmfieldnames{1}).Snew.position.xstep;

%defining some more useful variables
avgspec = mean(stxmdata.(stxmfieldnames{1}).Snew.spectr,3);
spec288 = stxmdata.(stxmfieldnames{1}).Snew.spectr(:,:,3);
STXMbinmap = stxmdata.(stxmfieldnames{1}).Snew.binmap;
STXMlabelmat = stxmdata.(stxmfieldnames{1}).Snew.LabelMat;

STXMimage = mat2gray(avgspec);
STXMimage = imadjust(STXMimage);

STXMrows = size(STXMimage,1);
STXMcols = size(STXMimage,2);



%fixing duplicate particles found in adjacent fields, This is a good thing
%to work on, but probably low priority.  Most FOV's are within a single
%field, those that aren't can just be excluded without the loss of too much
%data.
% % % % % % totnumfields = semdata.stub1.fieldposinfo.fieldnum;
% % % % % %
% % % % % % for jj = 1:numyfields;
% % % % % %     colnums(jj,:) = ((numxfields*(jj-1)+1):numxfields*jj); % this makes a vector [1 2 3] or [4 5 6]
% % % % % % end
% % % % % %
% % % % % % if size(colnums,2) > 1
% % % % % %     workingimgs1 = cat(2,semdata.stub1.imgdata(:,:,colnums(1,1)),semdata.stub1.imgdata(:,:,colnums(1,2)));
% % % % % %     compositefixname = 'PICK TWO POINTS REPRESENTING IDENTICAL PARTICLES (horizontal); IF NONE PRESENT, CLICK TWICE IN SAME SPOT';
% % % % % %     compositefix_fig1 = figure('name',compositefixname,...
% % % % % %         'Units','normalized',...
% % % % % %         'NumberTitle','off',...
% % % % % %         'Position',[0 0 1 1]);
% % % % % %     imshow(workingimgs1);
% % % % % %     [colfixpts_horiz, rowfixpts_horiz] = getpts(compositefix_fig1);
% % % % % %     close(compositefix_fig1);
% % % % % %     clear compositefix_fig1;
% % % % % %     dist2dup_horiz = abs(colfixpts_horiz(1) - colfixpts_horiz(2));
% % % % % % end
% % % % % %
% % % % % % if size(colnums,1) > 1
% % % % % %     workingimgs2 = cat(1,semdata.stub1.imgdata(:,:,colnums(1,1)),semdata.stub1.imgdata(:,:,colnums(2,1)));
% % % % % %     compositefixname2 = 'PICK TWO POINTS REPRESENTING IDENTICAL PARTICLES (vertical); IF NONE PRESENT, CLICK TWICE IN SAME SPOT';
% % % % % %     compositefix_fig2 = figure('name',compositefixname2,...
% % % % % %         'Units','normalized',...
% % % % % %         'NumberTitle','off',...
% % % % % %         'Position',[0 0 1 1]);
% % % % % %     imshow(workingimgs2);
% % % % % %     [colfixpts_vert, rowfixpts_vert] = getpts(compositefix_fig2);
% % % % % %     close(compositefix_fig2);
% % % % % %     clear compositefix_fig2;
% % % % % %     dist2dup_vert = abs(rowfixpts_vert(1) - rowfixpts_vert(2));
% % % % % % end
% % % % % %
% % % % % %
% % % % % % newcompositepic = []; %needed for cat command to work for the first matrix
% % % % % % for j = 1:numyfields
% % % % % %     rowimgdata = [];
% % % % % %     rownumvector = ((numxfields*(j-1)+1):numxfields*j); %this will create a vector of numbers corresponding to a single row (1 2 3) or (4 5 6) etc
% % % % % %     for q = 1:length(rownumvector)
% % % % % %         if q == 1
% % % % % %             rowimgdata = cat(2,rowimgdata,semdata.stub1.imgdata(:,:,rownumvector(q)));
% % % % % %         else
% % % % % %             rowimgdata = cat(2,rowimgdata,semdata.stub1.imgdata(:,dist2dup_horiz:end,rownumvector(q)));
% % % % % %         end
% % % % % %     end
% % % % % %     if j == 1
% % % % % %         newcompositepic = cat(1,newcompositepic,rowimgdata);
% % % % % %     else
% % % % % %         newcompositepic = cat(1,newcompositepic,rowimgdata(dist2dup_vert:end,:));
% % % % % %     end
% % % % % % end
% % % % % %
% % % % % % SEMimageraw_new = mat2gray(newcompositepic);
% % % % % % SEMimage = SEMimageraw_new;
% % % % % % SEMimage = imadjust(SEMimage,[0;1],[0;1],1);




% % % % % newcompositepic = semdata.stub1.imgdata(:,:,1); %base for which to build on
% % % % % for q = 1:max(max(colnums));
% % % % %
% % % % %     if any((q+1) == colnums(2:end,1)) == 1 %checks to see if a vertical match is needed
% % % % %         workingimgs = cat(1,semdata.stub1.imgdata(:,:,(q+1-numxfields)),semdata.stub1.imgdata(:,:,(q+1)));
% % % % %
% % % % %     else %else we do a horizontal matching
% % % % %         workingimgs = cat(2,semdata.stub1.imgdata(:,:,q),semdata.stub1.imgdata(:,:,q+1));
% % % % %         compositefix_fig = figure;
% % % % %         imshow(workingimgs);
% % % % %         [colfixpts, rowfixpts] = getpts(compositefix_fig);
% % % % %         dist2dup = abs(colfixpts(1) - colfixpts(2));
% % % % %         fixedimg = semdata.stub1.imgdata(:,dist2dup:end,(q+1));
% % % % %         newcompositepic = cat(2,newcompositepic,fixedimg);
% % % % %     end
% % % % %
% % % % % end

% SEMimage(SEMimage < 0.01) = 0.3;
%some spots are so dark that they are given a '0' intensity even though
%they are part of the particle.  This throws off automated particle
%detection efforts.  Replacing these dark spots with a number more
%"particle-like" hopefully will help this



stxmresizeratio = stxmpixsize./sempixsize;


%resizing STXM image based on pixel size
STXMresize = imresize(STXMimage,stxmresizeratio);
STXMresizebin = imresize(STXMbinmap,stxmresizeratio);
STXMlabelmatresize = imresize(STXMlabelmat,stxmresizeratio);


cd(stxmdir{1});
if exist([stxmfolder 'tformvars.mat'],'file') == 2
    usetformans = inputdlg('use saved tform and correlation data?','Load dialogue',1,{'yes'});
else
    usetformans = 'no';
end

if strcmp(usetformans,'yes')==1
    cd(stxmdir{1});
    load([stxmfolder 'tformvars.mat']);
    corrtest = 'yes';
else
    %2D correlation, finding where overlap is strongest
    corrimg = normxcorr2(STXMresize,SEMimage);
    [ypeak, xpeak] = find(corrimg==max(corrimg(:)));
    yoffset = ypeak- size(STXMresize,1);
    xoffset = xpeak - size(STXMresize,1);
    
    %position of correlation rectangle
    xstart = xoffset+1;
    xend = xstart + size(STXMresize,2);
    ystart = yoffset+1;
    yend = ystart + size(STXMresize,1);
    
    matchfig = figure;
    matchax = axes;
    imshow(SEMimage,'Parent',matchax);
    imrect(matchax, [xoffset+1, yoffset+1, size(STXMresize,2), size(STXMresize,1)]);
    movegui(matchfig,'northwest');
    matchfig.Name = 'Region found by correlation';
    matchfig.NumberTitle = 'off';
    
    
    % fstxmnorm = figure;
    % imshow(STXMimage);
    
    fcor = figure;
    surf(corrimg);
    shading flat
    movegui(fcor,'southeast');
    fcor.Name = 'Correlation Output';
    fcor.NumberTitle = 'off';
    
    buffersize = 10; %number of extra pixels to take to allow aspect ratio resizing
    
    
    %Useful definitions
    clipcolstart = xstart-buffersize;
    clipcolend = xend+buffersize;
    cliprowstart = ystart-buffersize;
    cliprowend = yend+buffersize;
    
    %taking the identified section and adding some padding
    if cliprowstart >= 1 && clipcolstart >= 1
        if cliprowend < size(SEMimage,1) && clipcolend < size(SEMimage,2);
            SEMclip = SEMimage((cliprowstart):(cliprowend),(clipcolstart):(clipcolend));
            
            fstxm = figure;
            fstxm.Units = 'normalized';
            fstxm.OuterPosition = [0.5 0 0.5 1];
            imshow(STXMresize);
            
            fsem = figure;
            fsem.Units = 'normalized';
            fsem.OuterPosition = [0 0 0.5 1];
            imshow(SEMclip);
            
            corrtest = inputdlg('Did the correlation work? (yes/no)','correlation test',1,{'yes'});
            close(fsem); clear fsem;
            close(fstxm); clear fstxm;
        else
            corrtest = inputdlg('Did the correlation work? (yes/no)','correlation test',1,{'yes'});
        end
    else
        corrtest = inputdlg('Did the correlation work? (yes/no)','correlation test',1,{'yes'});
    end
    
    %TEMPORARY FOR MSER OPTIMIZATION
    if cliprowstart < 1
        cliprowstart = 1;
    end
    
    if clipcolstart < 1;
        clipcolstart = 1;
    end
    
    if cliprowend > size(SEMimage,1)
        cliprowend = size(SEMimage,1);
    end
    
    if clipcolend > size(SEMimage,2)
        clipcolend = size(SEMimage,2);
    end

end

SEMclip = SEMimage((cliprowstart):(cliprowend),(clipcolstart):(clipcolend));
% corrtest = 'yes';

if strcmp(corrtest,'yes') == 1
else
    fsemimg = figure;
    fsemimg.Name = 'CHOOSE THE RECTANGLE THAT INCLUDES ALL OF THE STXM IMAGE, plus a small buffer';
    fsemimg.NumberTitle = 'off';
    fsemimg.Units = 'normalized';
    fsemimg.OuterPosition = [0 0 0.5 1];
    imshow(SEMimage);
    
    fstxm = figure;
    fstxm.Units = 'normalized';
    fstxm.OuterPosition = [0.5 0 0.5 1];
    imshow(STXMresize);
    
    manualrect = getrect(fsemimg); %manual rect is [xmin ymin width height]
    
    cliprowstart = round(manualrect(2));
    cliprowend = round(cliprowstart + manualrect(4));
    clipcolstart = round(manualrect(1));
    clipcolend = round(clipcolstart + manualrect(3));
    SEMclip = SEMimage((cliprowstart):(cliprowend),(clipcolstart):(clipcolend));
    
    close(fsemimg); clear fsemimg;
    close(fstxm); clear fstxm;
    
end

%taking out the data relevant to the clipped SEM image
cliprowidxs = find(semdata.stub1.semdata(:,6)>=cliprowstart  &  semdata.stub1.semdata(:,6)<=cliprowend);
clipcolidxs = find(semdata.stub1.semdata(:,7)>=clipcolstart  &  semdata.stub1.semdata(:,7)<=clipcolend);
clipidx = intersect(clipcolidxs,cliprowidxs);
clipsemdata = semdata.stub1.semdata(clipidx,:);

clipsemdata(:,6) = round(clipsemdata(:,6) - cliprowstart + 1); %the +1 is to ensure no 0's
clipsemdata(:,7) = round(clipsemdata(:,7) - clipcolstart + 1);

%finding and listing which elements we have data for
elementlist_start = find(strcmp(semdata.stub1.colheaders,'Orient'))+1;
elementlist_end = find(strcmp(semdata.stub1.colheaders,'AvgVideo'))-2;
elementstrlist = semdata.stub1.colheaders((elementlist_start):(elementlist_end));
numelements = length(elementstrlist);

for i = 1:numelements
    elementstrlist{i} = strtrim(elementstrlist{i}(1:2));
end

% % % % % %changing the pixel positions in "clipsemdata" into row/column positions
% % % % % clipcoldata = clipsemdata(:,6) - clipcolstart + 1; %columns start at 1 not zero, thus the +1
% % % % % cliprowdata = (size(SEMimage,1) - clipsemdata(:,7) + 1) - ystart + buffersize + 1;
% % % % %
% % % % % clipsemdata(:,6) = cliprowdata; %NOTE in colheaders, column 6 is x_cent and 7 is y_cent
% % % % % clipsemdata(:,7) = clipcoldata; %after this transformation, column 6 is columns and 7 is rows
% % % % %
% % % % % clipcolheaders = semdata.stub1.colheaders;
% % % % % clipcolheaders{6} = 'row_cent';
% % % % % clipcolheaders{7} = 'col_cent';


% % GrayImage=SEMclip; %mat2gray(imagebuffer); %% Turn into a greyscale with vals [0 1]
% % GrayImage=imadjust(GrayImage);%,[0 1],[0 1],15); %% increase contrast
% % Thresh=graythresh(GrayImage); %% Otsu thresholding
% % SEMclipbin=im2bw(GrayImage,Thresh); %% Give binary image
% % smallcorr = normxcorr2(STXMresize,SEMclip);
% % [ypeaksm, xpeaksm] = find(smallcorr==max(smallcorr(:)));
% % yoffsetsm = ypeaksm - size(STXMresize,1);
% % xoffsetsm = xpeaksm - size(STXMresize,1);
% % basecorrmax = max(max(smallcorr(:)));

% % fcor2 = figure;
% % surf(smallcorr);
% % shading flat
% % movegui(fcor2,'southwest');

% % matchfigsm = figure;
% % matchaxsm = axes;
% % imshow(SEMclip,'Parent',matchaxsm);
% % imrect(matchaxsm, [xoffsetsm+1, yoffsetsm+1, size(STXMresize,2), size(STXMresize,1)]);
% % movegui(matchfigsm,'northwest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(usetformans,'yes')==1;
    tforminv = invert(tform);
    SEMoutputview = imref2d(size(SEMclip));
    STXMxform = imwarp(STXMresize,tform,'OutputView',SEMoutputview);
    STXMbinxform = imwarp(STXMresizebin,tform,'OutputView',SEMoutputview);
    STXMlabelxform = imwarp(STXMlabelmatresize,tform,'OutputView',SEMoutputview);
        
    manualflag = 1;

    
    
else
    % MSER feature detecion
    regionsSEM = detectMSERFeatures(SEMclip,...
        'ThresholdDelta',3,...
        'RegionAreaRange', [1 14000],...
        'MaxAreaVariation',0.4);
    regionsSTXM = detectMSERFeatures(STXMresize,...
        'ThresholdDelta',3,...
        'RegionAreaRange', [1 14000],...
        'MaxAreaVariation',0.4);
    
    % Retrieving MSER regions and putting them in variables matlab uses
    [featuresSEM, validptsSEM] = extractFeatures(SEMclip,regionsSEM);
    [featuresSTXM, validptsSTXM] = extractFeatures(STXMresize,regionsSTXM);
    
    semptsnum = length(validptsSEM);
    stxmptsnum = length(validptsSTXM);
    
    fsem = figure;          imshow(SEMclip);    hold on;
    plot(regionsSEM)%, 'showPixelList', true, 'showEllipses', false);
    movegui(fsem,'southwest');
    
    fvalidsem = figure;     imshow(SEMclip);    hold on;
    plot(validptsSEM,'showOrientation',true);
    title('SEM');
    movegui(fvalidsem,'southeast');
    fvalidsem.Name = 'Points Found';
    fvalidsem.NumberTitle = 'off';
    
    fstxm = figure;     imshow(STXMresize);     hold on;
    plot(regionsSTXM)%, 'showPixelList', true, 'showEllipses', false);
    movegui(fstxm,'northwest');
    
    fvalidstxm = figure;    imshow(STXMresize); hold on;
    plot(validptsSTXM,'showOrientation',true);
    title('STXM');
    movegui(fvalidstxm,'northeast');
    
    % First pass at matching the two extracted MSER features between images
    indexpairs = matchFeatures(featuresSEM, featuresSTXM,...
        'MatchThreshold', 90,...
        'MaxRatio',0.8,...
        'Unique', false);
    
    % Retrieving list of only the matched features
    matchedSEM = validptsSEM(indexpairs(:,1));
    matchedSTXM = validptsSTXM(indexpairs(:,2));
    
    matchnum = length(indexpairs);
    
    fmatch = figure;
    showMatchedFeatures(SEMclip, STXMresize,matchedSEM,matchedSTXM);
    movegui(fmatch,'southwest');
    fmatch.Name = 'Matched Features';
    fmatch.NumberTitle = 'off';
    
    %Defining the coordinate system imwarp will use
    SEMoutputview = imref2d(size(SEMclip));
    
    if length(matchedSEM) >= 4 || length(matchedSTXM) >= 4;
        % Determining/excluding outliers AND calculating the transform matrix
        % relating the two images  SEMimage = tform(STXMimage)
        [tform, inlierSTXM, inlierSEM] = estimateGeometricTransform(...
            matchedSTXM,matchedSEM, 'projective',...
            'MaxNumTrials',50000,...
            'Confidence', 99,...
            'MaxDistance',16);
        
        inliernum = length(inlierSEM);
        
        finliers = figure;
        showMatchedFeatures(SEMclip, STXMresize, inlierSEM, inlierSTXM);
        movegui(finliers,'northwest');
        finliers.Name = 'Inlier Matches';
        finliers.NumberTitle = 'off';
        
        Tinv = tform.invert.T; %this is nifty but isn't a "transform variable" matlab can use on images, only good for retrieving scale and angle.
        
        %Retrieving the inverse transform STXMimage = tforminv(SEMimage);
        tforminv = invert(tform);
        
        
        %some possibly useful varibles
        ss = Tinv(2,1);
        sc = Tinv(1,1);
        scale_rec = sqrt(ss*ss + sc*sc);
        theta_rec = atan2(ss,sc)*180/pi;
        
        %Defining the coordinate system imwarp will use
        SEMoutputview = imref2d(size(SEMclip));
        
        %Applying the tform matrix, (rotation and scale change)
        STXMxform = imwarp(STXMresize,tform,'OutputView',SEMoutputview);
        
        pairf = figure;
        imshowpair(SEMclip,STXMxform);
        movegui(pairf,'north');
        pairf.Name = 'Overlaid Images';
        pairf.NumberTitle = 'off';
        
    end
    
    
    qualitycheck = inputdlg('Did the automatic method work (and continue) or no (and match manually)? [yes/no]','auto/manual check',1,{'yes'});
    
    if strcmp(qualitycheck,'yes') == 1
        
        % Transforming bin and label STXM mats
        STXMbinxform = imwarp(STXMresizebin,tform,'OutputView',SEMoutputview);
        STXMlabelxform = imwarp(STXMlabelmatresize,tform,'OutputView',SEMoutputview);
        
        manualflag = 0;
        
    else
        clear STXMxform tform tforminv
        
        figname_help = 'CHOOSE n PTS FROM IMAGE THEN PRESS ENTER';
        manualmatchingfigSEM = figure('name',figname_help,...
            'NumberTitle','off');
        manualmatchingfigSEM.Units = 'normalized';
        manualmatchingfigSEM.OuterPosition = [0 0 0.5 1];
        
        imshow(SEMclip);
        title('SEM clip');
        
        manualmatchingfigSTXM = figure('name',figname_help,...
            'NumberTitle','off');
        manualmatchingfigSTXM.Units = 'normalized';
        manualmatchingfigSTXM.OuterPosition= [0.5 0 0.5 1];
        imshow(STXMresize);
        title('STXM resize');
        
        %     numpts2select = inputdlg('How many points to select? (number only)','number of points selection',1,{'5'});
        
        %     helpprompt = sprintf('%s','Choose',numpts2select{1},'points from the left picture, then repeat for the right picture');
        %     helpdlg(helpprompt);
        
        help_hdl = helpdlg('Choose pts with L-click, end with double-click, R-click, or enter.  Delete/backspace to go back');
        movegui(help_hdl,'northeast')
        
        [SEMx_selected, SEMy_selected] = getpts(manualmatchingfigSEM);
        figure(manualmatchingfigSEM);
        hold on
        plot(SEMx_selected,SEMy_selected,'r+');
        
        %making number list so that selected points are numbered in order
        numlabels = cellstr(num2str((1:length(SEMx_selected))'));
        for kk = 1:length(SEMx_selected);
            text(SEMx_selected(kk),SEMy_selected(kk),numlabels{kk},...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','right',...
                'Color','r');
        end
        if exist('help_hdl','var') ;
            close(help_hdl);
        end
        
        help_hdl2 = helpdlg('Choose pts with L-click, end with double-click, R-click, or enter.  Delete/backspace to go back');
        movegui(help_hdl2,'northwest')
        
        [STXMx_selected, STXMy_selected] = getpts(manualmatchingfigSTXM);
        figure(manualmatchingfigSTXM);
        
        if exist('help_hdl2','var') ;
            close(help_hdl2);
        end
        %     close(manualmatchingfig); clear manualmatchingfig;
        %     lcombo = length(SEMx_selected_thenSTXMx);
        %     SEMx_selected = SEMx_selected_thenSTXMx(1:(lcombo/2));
        %     SEMy_selected = SEMy_selected_thenSTXMy(1:(lcombo/2));
        %     STXMx_selected = SEMx_selected_thenSTXMx((lcombo/2 + 1):lcombo);
        %     STXMy_selected = SEMy_selected_thenSTXMy((lcombo/2 + 1):lcombo);
        
        movingPoints = cat(2,STXMx_selected,STXMy_selected);
        fixedPoints = cat(2,SEMx_selected,SEMy_selected);
        
        % %     cpselect(STXMresize,SEMclip); %moving, fixed  This is so useful but it doesn't seem to want to work inside a function
        tform = fitgeotrans(movingPoints,fixedPoints,'projective');
        tforminv = invert(tform);
        
        STXMxform = imwarp(STXMresize,tform,'OutputView',SEMoutputview);
        STXMbinxform = imwarp(STXMresizebin,tform,'OutputView',SEMoutputview);
        STXMlabelxform = imwarp(STXMlabelmatresize,tform,'OutputView',SEMoutputview);
        
        pairf2 = figure;
        imshowpair(SEMclip,STXMxform);
        movegui(pairf2,'northeast');
        
        manualflag = 1;
        close(manualmatchingfigSEM);
        close(manualmatchingfigSTXM);
        
    end
end
% Size definitions
[xformrow,xformcol] = size(STXMlabelxform);
[labelrow,labelcol] = size(STXMlabelmat);

% Defining another coordinate system
resizeoutputview = imref2d(size(STXMlabelmatresize)); %I'm still not quite sure what this does

% This is part of being able to move completely from STXM to SEM
% coordinates
stxmresizeratio_reverse = 1./stxmresizeratio;

% Some preallocation
zerosmap = zeros(xformrow,xformcol);
EDXmap = zeros(labelrow,labelcol,numelements);
nostxm_parts = zeros(10,2);

% This takes each pixel of a binary map, and determines the distance to the
% nearest non-zero (Dist2partsmap) AND gives the linear index of which
% non-zero pixel is the closest (lindexofnearestpartmap).  Both are
% returned as matricies of the same size as the input matrix.
[Dist2partsmap,lindexofnearestpartmap] = bwdist(STXMbinmap);

if any(exist('sillystring','file'))
    hwait = waitbar(0,sillystring);
else
    hwait = waitbar(0,'plz w8');
end

nostxm_idx = 1;

%loops over number of elements EDX found
for k = 1:numelements
    %each loop makes a map of zeros with a single 1 at the center of the
    %"i"th particle (determined by CCSEM).  This map is in the SEM
    %coordinate system.  This same map is transformed and resized into the
    %STXM coordinate system.  It then checks the position and finds the
    %particle number (as found in "STXMlabelmat").  This then proceeds to
    %make every position in STXMlabelmat (with that particle number) equal
    %to the relative % of a given element (the "k"th element).  This gives
    %a map where each particle is assumed to have a homogeneous spread of
    %the "k"th element.
    for i = 1:size(clipsemdata,1) %loops over number of particles determined by CCSEM
        partlabel = zerosmap; %rezeros this blank map
        partlabel(clipsemdata(i,6),clipsemdata(i,7)) = 1; %makes the center of the particle (found by CCSEM)
        partinverted = imwarp(partlabel,tforminv,'OutputView',resizeoutputview);
        partsize_reverse = imresize(partinverted,stxmresizeratio_reverse);
        peakval = max(max(partsize_reverse));  %the transforms turn the single pixel "1" value into a small blob,
        %the peak is assumed to be the correct position, this is close enough
        [peak_row, peak_col] = find(partsize_reverse == peakval,1); %this finds the row/col position of the above maximum. This could probably be
        %incorporated into the output of the "max" function but then it would give
        %linear index which I then would need to split.  It's good like this I think.
        
        if peak_row > size(STXMlabelmat,1) || peak_col > size(STXMlabelmat,2); %if the SEM particle "i" is outside the bounds of the STXM image
            nostxm_parts(nostxm_idx,1) = peak_row;
            nostxm_parts(nostxm_idx,2) = peak_col;
            nostxm_idx = nostxm_idx + 1;
            continue
        end
        
        labelmat_partnum = STXMlabelmat(peak_row,peak_col); %finding the particle label number at the determined position
        
        if labelmat_partnum ~= 0 %if a particle position number is found
            part_lindex2d = find(STXMlabelmat == labelmat_partnum); %row/column format gives boxes instead of shaped particles
            part_lindex3d = labelrow.*labelcol.*(k-1) + part_lindex2d;  %this adds the number of cells in each matrix to allow the
            %2d lindex to be used in the 3d EDXmap matrix
            EDXmap(part_lindex3d) = clipsemdata(i,(elementlist_start+k-1)); %applies the "k"th element data to the position(s) found at part_lindex
            
        elseif labelmat_partnum == 0  %if the SEM particle "i" wasn't picked up with STXM or the particle finding routine
            if manualflag == 0;
                if Dist2partsmap(peak_row,peak_col) < 10; %if nearest particle is within 10 pixels, apply that data to the close particle instead of dumping it
                    closepartlindex = lindexofnearestpartmap(peak_row,peak_col);
                    labelmat_partnum = STXMlabelmat(closepartlindex);
                    part_lindex2d = find(STXMlabelmat == labelmat_partnum); %row/column format gives boxes instead of shaped particles
                    part_lindex3d = labelrow.*labelcol.*(k-1) + part_lindex2d;  %this adds the number of cells in each matrix to allow the
                    %2d lindex to be used in the 3d EDXmap matrix
                    EDXmap(part_lindex3d) = clipsemdata(i,(elementlist_start+k-1));
                    
                else %if no near particles, SEM particle "i" wasn't picked up with STXM
                    nostxm_parts(nostxm_idx,1) = peak_row;
                    nostxm_parts(nostxm_idx,2) = peak_col;
                    nostxm_idx = nostxm_idx + 1;
                end
            elseif manualflag == 1;
                if Dist2partsmap(peak_row,peak_col) < 15; %if images are manually aligned, extra buffer is given for loss of precision selecting two points by hand
                    closepartlindex = lindexofnearestpartmap(peak_row,peak_col);
                    labelmat_partnum = STXMlabelmat(closepartlindex);
                    part_lindex2d = find(STXMlabelmat == labelmat_partnum); %row/column format gives boxes instead of shaped particles
                    part_lindex3d = labelrow.*labelcol.*(k-1) + part_lindex2d;  %this adds the number of cells in each matrix to allow the
                    %2d lindex to be used in the 3d EDXmap matrix
                    EDXmap(part_lindex3d) = clipsemdata(i,(elementlist_start+k-1));
                    
                else %if no near particles, SEM particle "i" wasn't picked up with STXM
                    nostxm_parts(nostxm_idx,1) = peak_row;
                    nostxm_parts(nostxm_idx,2) = peak_col;
                    nostxm_idx = nostxm_idx + 1;
                end
            end
        end
    end
    waitbar(k./numelements,hwait);
end
close(hwait);



Cmassmap = stxmdata.(stxmfieldnames{1}).Snew.elemap.C .* stxmdata.(stxmfieldnames{1}).Snew.binmap;
Cerrmassmap = stxmdata.(stxmfieldnames{1}).Snew.errmap.C .* stxmdata.(stxmfieldnames{1}).Snew.binmap;
Nmassmap = stxmdata.(stxmfieldnames{1}).Snew.elemap.N .* stxmdata.(stxmfieldnames{1}).Snew.binmap;
Nerrmassmap = stxmdata.(stxmfieldnames{1}).Snew.errmap.N .* stxmdata.(stxmfieldnames{1}).Snew.binmap;
Omassmap = stxmdata.(stxmfieldnames{1}).Snew.elemap.O .* stxmdata.(stxmfieldnames{1}).Snew.binmap;
Oerrmassmap = stxmdata.(stxmfieldnames{1}).Snew.errmap.O .* stxmdata.(stxmfieldnames{1}).Snew.binmap;

[EDXrow, EDXcol, EDXele] = size(EDXmap);

oxystr = strfind(elementstrlist,'O');
for q = 1:length(oxystr);
    oxytest = any(oxystr{q});
    if oxytest == 1;
        oxyidx = q;
        break
    end
end

if oxyidx == 2 %Oxygen is either at position 2 or 3 depending on if Nitrogen was taken, Nitrogen will always be at position 2
    EDXmap_new = zeros(EDXrow,EDXcol,EDXele+1);
    EDXmap_new(:,:,1) = Cmassmap;
    EDXmap_new(:,:,2) = Nmassmap;
    EDXmap_new(:,:,3) = Omassmap;
    
    elementstrlist_new = cell(1,EDXele+1);
    elementstrlist_new{1} = elementstrlist{1};
    elementstrlist_new{2} = 'N';
    elementstrlist_new{3} = 'O';
    
    for x = 4:(EDXele+1)
        EDXmap_new(:,:,x) = EDXmap(:,:,x-1);
        elementstrlist_new{x} = elementstrlist{x-1};
    end
    EDXmap = EDXmap_new;
    elementstrlist = elementstrlist_new;
elseif oxyidx == 3
    EDXmap(:,:,1) = Cmassmap;
    EDXmap(:,:,2) = Nmassmap;
    EDXmap(:,:,3) = Omassmap;
end

EDXerrmap = zeros(size(EDXmap));
EDXerrmap(:,:,1) = Cerrmassmap;
EDXerrmap(:,:,2) = Nerrmassmap;
EDXerrmap(:,:,3) = Oerrmassmap;

for gg = 4:size(EDXmap,3)
    EDXerrmap(:,:,gg) = EDXmap(:,:,gg) .* 0.05; %error estimation of about 5%
end

Alstr = strfind(elementstrlist,'Al');
for qq = 1:length(Alstr);
    Altest = any(Alstr{qq});
    if Altest == 1;
        Alidx = qq;
        EDXmap(:,:,qq) = [];
        break
    end
end

elementstrlist_temp = elementstrlist;
for h = qq:(length(elementstrlist)-1);
    elementstrlist_temp{h} = elementstrlist{h+1};
end
elementstrlist = elementstrlist_temp;

new_numelements = size(EDXmap,3);
%%%%%%%%%Plots all element data in a multi-plot image
EDXfigure = figure;
set(gcf,'Units','normalized','Position',[0.005 0.05 0.99 0.86]);
subploth = tight_subplot(3,5,[.04 .01],[.03 .03],[.03 .03]);
for j = 1:new_numelements;
%     subplot(3,5,j);
    axes(subploth(j));
    imagesc(EDXmap(:,:,j));
    set(subploth,'xticklabel',[],'yticklabel',[]);
    title(elementstrlist{j});
%     axis square
    colorbar;
end
hold on
% subplot(3,5,15);
axes(subploth(15));
imagesc(STXMbinmap);
set(gca,'xticklabel',[],'yticklabel',[]);
title('STXM binmap');
% axis square
colorbar;

%If you're using saved tform data, no need to save it again.
if strcmp(usetformans,'yes') == 1
    savetformanswer = 'no';
else
    savetformanswer = inputdlg('Do you want to save correlation and tform data (so you don''t have to pick regions again)?','Save dialogue',1,{'yes'});
end

if strcmp(savetformanswer,'yes') == 1
    cd(stxmdir{1});
    tformvarsfilename = [stxmfolder 'tformvars'];
    save(tformvarsfilename,'cliprowstart','cliprowend','clipcolstart','clipcolend','tform');
end

inputanswer = inputdlg('Did it work? Proceed and save stuff?(yes or no)','Save dialogue',1,{'yes'});

if strcmp(inputanswer,'yes') == 1
    cd(stxmdir{1});
    figurenametest = sprintf('%s',stxmfolder,'_EDX*');
    figurename = sprintf('%s',stxmfolder,'_EDX');
    if isempty(dir(figurename))
        savefig(EDXfigure,figurename);
        saveas(EDXfigure,figurename,'png');
        saveas(EDXfigure,figurename,'eps');
        saveas(EDXfigure,figurename,'tif');
    end
    
    varfilename = sprintf('%s',stxmfolder,'_EDXvars');
    varfilenamecheck = sprintf('%s',varfilename,'.mat');
    if isempty(dir(varfilenamecheck))
        clear pairf finliers fmatch fvalidsem fvalidstxm fsem SEMoutputview
        clear fstxm fcor fstxmnorm fstxm EDXfigure hwait matchax matchfig
        clear manualma* subploth tform tforminv validpts*
        clear pairf2 ans help_* matched* regions* resizeoutputview
        save(stxmfolder);
    end
    
else
end


% % semclipf = figure;
% % imshow(SEMclip);
% % movegui(semclipf,'northwest');
% %
% % stxmxformf = figure;
% % imshow(STXMxform);
% % movegui(stxmxformf,'northeast');




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


end
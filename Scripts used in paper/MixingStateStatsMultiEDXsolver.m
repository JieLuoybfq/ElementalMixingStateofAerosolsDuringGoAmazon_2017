
%Determining mass fractions without


%this part was removed to allow it to be run much faster without picking
%data set each time
prompt = {'How many data sets to be analyzed? (number only)'};
numsetsSTR = inputdlg(prompt,'Input',1);
numsets = str2double(numsetsSTR);
numsets = 9;

%this part was removed to allow it to be run much faster without picking
%data set each time
datasettitle = cell(numsets,1);
filedirsSET = cell(numsets,1);
options.Resize = 'on';
options.WindowStyle = 'normal';
for i = 1:numsets
	prompt2 = sprintf('%s','What would you like to title data set #',numsets,'? (matlab variable rules)');
	datasettitle{i} = inputdlg(prompt2,'Data Set Name',1,{''},options);
	filedirsSET{i} = uipickfiles; %ui window for picking multiple directories
end

% load('EDXfiledirsSET_home_NoOxyremoved_bigZF2removed.mat');

% datasettitle = cell(9,1);
% % datasettitle{1}{1} = 'Atto14';
% datasettitle{1}{1} = 'Atto15s7';
% datasettitle{2}{1} = 'T313ds7';
% datasettitle{3}{1} = 'T313ds8';
% datasettitle{4}{1} = 'T313ns7';
% datasettitle{5}{1} = 'T314ds7';
% datasettitle{6}{1} = 'T314ds8';
% datasettitle{7}{1} = 'ZF2W20s7';
% datasettitle{8}{1} = 'ZF2W21s7';
% datasettitle{9}{1} = 'ZF2W21s8';
% 
% filedirsSET = cell(10,1);
% filedirsSET{1} = cell(1,1);
% filedirsSET{1}{1} = 'C:\Users\Emily\Google Drive\Projects\GoAmazon\All Amazon STXM Data\Atto Oct 14\150618 Beamtime\ShortStack532_150620015';
% filedirsSET{2} = cell(9,1);
% filedirsSET{2}{1} = 'C:\Users\Emily\Google Drive\Projects\GoAmazon\All Amazon STXM Data\Atto Oct 15 stg 7\150618 Beamtime\ShortStack532_150619075';
% filedirsSET{2}{2} = '
% filedirsSET{2}{3} = '
% filedirsSET{2}{4} = '
% filedirsSET{2}{5} = '
% filedirsSET{2}{6} = '
% filedirsSET{2}{7} = '
% filedirsSET{2}{8} = '
% filedirsSET{2}{9} = '



% clear options;


mfracelelist = cell(14,1);
mfracelelist{1} = 'C';
mfracelelist{2} = 'N';
mfracelelist{3} = 'O';
mfracelelist{4} = 'Na';
mfracelelist{5} = 'Mg';
mfracelelist{6} = 'P';
mfracelelist{7} = 'S';
mfracelelist{8} = 'Cl';
mfracelelist{9} = 'K';
mfracelelist{10} = 'Ca';
mfracelelist{11} = 'Mn';
mfracelelist{12} = 'Fe';
mfracelelist{13} = 'Ni';
mfracelelist{14} = 'Zn';

ulist = STXMEDXulist; %this is a 14x6 cell containing mass absorption coefficients for the 14 elements listed above




for ee = 1:numsets %looping through each directory SET
    clear filedirs MixingOverview %this is to stop a larger fileset (9 FOV's) to keep the last few FOV's while the new, smaller set, populates only 1-6 and leaves 7,8,9 alone
    filedirs = filedirsSET{ee};
    
    totalmfracmatrix = [];
    totalDimatrix = [];
    
    tempDa = zeros(length(filedirs),1);
    tempDy = tempDa;
    tempDb = tempDa;
    tempCHI = tempDa;
    tempnumparts = tempDa;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%Dealing with a specific data set%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if any(exist('sillystring','file'))
		hwait = waitbar(0,[sillystring ' (' num2str(ee) '/' num2str(numsets) ')']);
	else
		hwait = waitbar(0,'plz w8');
    end
    
    for jj = 1:length(filedirs); %looping through directories selected in the larger set
                
        [~,foldername,~] = fileparts(filedirs{jj});
        cd(filedirs{jj});
        matfilename = sprintf('%s',foldername,'.mat');
        tformfilename = [foldername, 'tformvars.mat'];
        clear elementstrlis* %this is to make sure whenever a new matfilename is loaded, a new elementstrlist and possibly a new elementstrlist_new is loaded fresh
        if ~isempty(dir(matfilename))
            load(matfilename);
            load(tformfilename);
            tforminv = invert(tform);
        else
            errmsg = sprintf('%s',matfilename,'doesn"t exist');
            continue
        end
        
        
        %this is needed or else this script will use too much memory
        clear avgspec buffersize clipcol* cliprow* clipidx clipsemdata closepartlindex coldata corrimg
        clear corrtest Dist2partsmap featuresSEM featuresSTXM fieldxsize fieldysize figname_help figurename figurenametest
        clear fixedPoints help_hdl help_hdl2 inlierS* labelcol* labelmat_partnum labelrow lindexofnearest*
        clear manualflag matchedS* movingPoints* nostxm* part_lindex* partinverted partlabel partsize_reverse peak_col peak_row
        clear peakval qualitycheck regionsS* resizeoutputview rowdata sc scale_rec SEMimage SEMimageraw spec288
        clear STXMbinmap STXMbinxform STXMcols STXMimage STXMlabelmat STXMlabelmatresize STXMlabelxform 
        clear STXMresize STXMresizebin STXMrows STXMx_selected STXMxform STXMy_selected theta_rec Tinv xfields zerosmap
        clear validptsS* xend xformcol xformrow xoffset xpeak xstart yend yfields yoffset ypeak ystart 
        clear elementl* field* fig* fixed* indexp* moving* nostxm*
        clear numx* numy* quality* sc* semp* semd* stxmp*
        
        EDXmap_old = EDXmap; %this saves the old EDXmap and allows me to replace EDXmap with the new EDXmap so the rest of the code still works
        [EDXrow, EDXcol, EDXele] = size(EDXmap_old);
        rgbhandle = stxmdata.(stxmfieldnames{1}).Snew.RGBCompMap;
        [crow,ccol,~] = size(rgbhandle);
        tempRGBmat = zeros(crow,ccol);
        tempRGBmat(rgbhandle(:,:,1)==255) = 1;
        tempRGBmat(rgbhandle(:,:,2)==170) = 2;
        tempRGBmat(rgbhandle(:,:,3)==255) = 3;
        
        
        for j = 1:length(mfracelelist);
            eleidx2 = find(strcmp(mfracelelist{j},elementstrlist)==1); %this finds the index in elementstrlist (and thus also in EDXmap) of the jth element in mfracelelist (the complete list)
            if ~isempty(eleidx2)
                EDXmap(:,:,j) = EDXmap_old(:,:,eleidx2(1));
            else
                EDXmap(:,:,j) = zeros(EDXrow,EDXcol);
            end
        end
        
        stxmfields = fieldnames(stxmdata);
        Cmassmap = stxmdata.(stxmfields{1}).Snew.elemap.C .* stxmdata.(stxmfields{1}).Snew.binmap;
        Nmassmap = stxmdata.(stxmfields{1}).Snew.elemap.N .* stxmdata.(stxmfields{1}).Snew.binmap;
        Omassmap = stxmdata.(stxmfields{1}).Snew.elemap.O .* stxmdata.(stxmfields{1}).Snew.binmap;
        
        EDXmap(:,:,1) = Cmassmap;
        EDXmap(:,:,2) = Nmassmap;
        EDXmap(:,:,3) = Omassmap;
        
        
        
        %         Crelpercentmap = EDXmap(:,:,1);
        %         Convratiomap = Cmassmap./Crelpercentmap; %making conversion ratio based on C.  Each pixel is it's own ratio
        %         nanmap = isnan(Convratiomap);
        %         Convratiomap(nanmap==1) = 0;
        %         Convratiomap(isinf(Convratiomap)==1) = 0;
        

        
        %labeling particles found by both SEM and STXM
        SEMlabelmap = bwlabel(sum(EDXmap_old,3)); %this sums up all elemental maps as determined by EDX and so has the added benefit of recognizing particles without any carbon
        numparticles = max(max(SEMlabelmap));
        
        
        avgCmap = zeros(EDXrow,EDXcol);
        avgNmap = zeros(EDXrow,EDXcol);
        avgOmap = zeros(EDXrow,EDXcol);
        avgODmap = zeros(EDXrow,EDXcol,size(stxmdata.(stxmfieldnames{1}).Snew.spectr,3));
        ODmap = stxmdata.(stxmfieldnames{1}).Snew.spectr;
        errODmap = stxmdata.(stxmfieldnames{1}).Snew.errOD;
        STXMenergy = stxmdata.(stxmfieldnames{1}).Snew.eVenergy;
        for j = 1:numparticles
            partmask = zeros(EDXrow,EDXcol);
            partmask(SEMlabelmap==j) = 1;
            currCpartmap = Cmassmap.*partmask;
            currNpartmap = Nmassmap.*partmask;
            currOpartmap = Omassmap.*partmask;
            
            currCavg = mean(mean(currCpartmap(currCpartmap>0))); %probably should have mean(mean(currCpartmap(currCpartmap>0))) otherwise it is finding the mean of a bunch of 0 pixels
            currNavg = mean(mean(currNpartmap(currNpartmap>0)));
            currOavg = mean(mean(currOpartmap(currOpartmap>0)));
            
            
            
            avgCmap(SEMlabelmap==j) = currCavg;
            avgNmap(SEMlabelmap==j) = currNavg;
            avgOmap(SEMlabelmap==j) = currOavg;
            
            
            for p = 1:length(STXMenergy);
                currODpart = ODmap(:,:,p) .* partmask;
                currODpart(currODpart <0) = 0;
                currerrODpart = errODmap(:,:,p) .* partmask;
                currODavg = mean(mean(currODpart(currODpart>0)));
                currerrODavg = std(reshape(currerrODpart(currerrODpart>0),[numel(currerrODpart(currerrODpart>0)) 1]))./sqrt(numel(currerrODpart(currerrODpart>0)));
                partlindex = find(SEMlabelmap == j); %linear index of where particle is in each 2D map
                avgODmap(partlindex + (p-1).*EDXrow.*EDXcol) = currODavg;
                errODmap(partlindex + (p-1).*EDXrow.*EDXcol) = currerrODavg;
            end
            
        end
        
        STXM_EDXmap = EDXmap; %this is a map of EDX elements with STXM maps for C/N/O with sub-particle resolution
        
        EDXmap(:,:,1) = avgCmap;
        EDXmap(:,:,2) = avgNmap;
        EDXmap(:,:,3) = avgOmap; %Now EDXmap has each particle with only a single value for either mass or relative mass fraction
        
        Elefracmap = zeros(size(EDXmap));
        for j = 1:numparticles
            partmask = zeros(EDXrow,EDXcol);
            partmask(SEMlabelmap==j) = 1;
            [partrow, partcol] = find(SEMlabelmap == j);
            
            
            cropvals = [min(partrow)   max(partrow)   min(partcol)   max(partcol)];
            
            for p = 1:length(mfracelelist)
                EDXmap(isnan(EDXmap)==1) = 0; %this replaces any NaN's found in EDXmap with 0's instead
                eletestmap = EDXmap(:,:,p).*partmask;
                eletest = sum(sum(eletestmap));
                if any(eletest)==1
                    eleflaglist(p) = 1;
                else
                    eleflaglist(p) = 0;
                end
                
                EDXmassvar(p) = EDXmap(partrow(1),partcol(1),p);
                
            end
            
            CO_EDXmass(j,jj,ee) = EDXmassvar(1)./EDXmassvar(3);
            CN_EDXmass(j,jj,ee) = EDXmassvar(1)./EDXmassvar(2);
            
            for q = 1:length(STXMenergy)
                if STXMenergy(q)>260 && STXMenergy(q)<282
                    ODlist(1) = avgODmap(partrow(1),partcol(1),q);
                    ERRODlist(1) = errODmap(partrow(1),partcol(1),q);
                elseif STXMenergy(q)>315 && STXMenergy(q)<325
                    ODlist(2) = avgODmap(partrow(1),partcol(1),q);
                    ERRODlist(2) = errODmap(partrow(1),partcol(1),q);
                elseif STXMenergy(q)>390 &&STXMenergy(q)<405
                    ODlist(3) = avgODmap(partrow(1),partcol(1),q);
                    ERRODlist(3) = errODmap(partrow(1),partcol(1),q);
                elseif STXMenergy(q)>425 &&STXMenergy(q)<435
                    ODlist(4) = avgODmap(partrow(1),partcol(1),q);
                    ERRODlist(4) = errODmap(partrow(1),partcol(1),q);
                elseif STXMenergy(q)>520 &&STXMenergy(q)<530
                    ODlist(5) = avgODmap(partrow(1),partcol(1),q);
                    ERRODlist(5) = errODmap(partrow(1),partcol(1),q);
                elseif STXMenergy(q)>545 &&STXMenergy(q)<555
                    ODlist(6) = avgODmap(partrow(1),partcol(1),q);
                    ERRODlist(6) = errODmap(partrow(1),partcol(1),q);
                end
               
                
            end
            ODlist(isnan(ODlist)==1) = 0;
            
            %%%%%%%%%%%%%
            %%%%Adding in error to max resulting chi, this is done by %%%%%
            %%%%%%%%%%
            %
            %for t = 1:6  %this switch case set will increase the carbon edge jump while decreaseing all others, effectively boosting amount of carbon
            %    switch t
            %        case 2 || 4 || 6
            %            ODlist(t) = ODlist(t) - ERRODlist(t);
            %       case 1 || 3 || 5
            %            ODlist(t) = ODlist(t) - ERRODlist(t);
            %    end
            %end
            %
            %EDXmassvar = EDXmassvar.*0.99;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%SOLVER STUFF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            startingfracs = zeros(1,14);
            startingfracs(1,1) = 0.4;%EDXmassvar(1)./sum(EDXmassvar(1:3)); %0.4;
            startingfracs(1,2) = 0.2;%EDXmassvar(2)./sum(EDXmassvar(1:3)); % 0.2;
            startingfracs(1,3) = 0.3;%EDXmassvar(3)./sum(EDXmassvar(1:3)); %0.3;
            
            startingfracs(1,4:14) = 0.05;
            startingfracs(1,15) = 0.00002; %0.000036
            
            fraclb = zeros(1,14);
            fracub = ones(1,14);
            %this should be rho*t = (0.5 - 2.5) * (0.000001 - 0.0005) where "-" represent ranges for rho and t respectively
            fraclb(1,15) = 0.5.*0.000001;  %used to be 0
            fracub(1,15) = 2.5.*0.0010;  %used to be 100
            
            dummyf = @(f) ODsystem(f,ODlist,EDXmassvar,eleflaglist); %using dummy function to allow ODsystem to be solved for f but take some parameters
            
            options = optimoptions('lsqnonlin','Display','none','TolFun',1e-9,'TolX',1e-9,'DiffMaxChange',0.1);
            [fsolutions,fval] = lsqnonlin(dummyf,startingfracs,fraclb,fracub,options);
            eleidx_zeros = find(eleflaglist==0);
            fsolutions(eleidx_zeros) = 0;
            
            savestartfracC(j,jj,ee) = startingfracs(1);
            savestartfracN(j,jj,ee) = startingfracs(2);
            savestartfracO(j,jj,ee) = startingfracs(3);
            savefsolutionsC(j,jj,ee) = fsolutions(1);
            savefsolutionsN(j,jj,ee) = fsolutions(2);
            savefsolutionsO(j,jj,ee) = fsolutions(3);
            savefsolutionsCO(j,jj,ee) = fsolutions(1)./fsolutions(3);
            savefsolutionsCN(j,jj,ee) = fsolutions(1)./fsolutions(2);
            %             if iscolumn(fsolutions)==1
            %                 fmatrix = fsolutions(1:14);
            %             else
            %                 fmatrix = fsolutions(1:14)';
            %             end
            
            %             if iscolumn(ODlist)==1
            %                 ansmatrix = ODlist;
            %             else
            %                 ansmatrix = ODlist';
            %             end
            %             ERRansmatrix = ERRODlist'; %when a vector is formed by var(1), var(2)... it is a row vector, this needs to be a column vec.
            %
            %
            %             ulist = STXMEDXulist; %this is a 14x6 matrix
            %             R = fsolutions(15);
            %             coeffmatrix = cell2mat(ulist') .*R; %This makes a 6x14 coefficient matrix
            %             ERRcoeffmatrix = coeffmatrix .* 0.05;
            %             coeffmatrix(7,:) = zeros(1,14);
            %             coeffmatrix(8,:) = zeros(1,14);
            %             %this mess is adding 2 equations to our system of unkn's and determines how
            %             %to get those two equations using C,N, and O data.  Either 2 useful
            %             %equations come out, or 2 trivial equations come out depending on if C,N,
            %             %and O data exists.
            %             eleidx2 = find(eleflaglist==1); %retrieves linear index of all non-zero elements
            %             if ~isempty(eleidx2(eleidx2 == 1)); %if eleidx(lindex of nonzero elements) is 1, C is present
            %                 if ~isempty(eleidx2(eleidx2 ==2)); %if eleidx(lindex of nonzero elements) is 2, N is present
            %                     %systemout(7) = -f(1) + f(2).*(mass(1)./mass(2));
            %                     coeffmatrix(7,2) = EDXmassvar(1)./EDXmassvar(2);
            %                     ERRcoeffmatrix(7,2) = sqrt(((EDXmassvar(1).*0.05)./EDXmassvar(2)).^2 + (EDXmassvar(1).*EDXmassvar(2).*0.05./EDXmassvar(2).^2).^2);
            %                     ansmatrix(7) = fsolutions(1);
            %                     if ~isempty(eleidx2(eleidx2 == 3)); %if eleidx(lindex of nonzero elements) is 3, O is present
            %                         %systemout(8) = -f(2) + f(3) .*(mass(2)./mass(3));
            %                         coeffmatrix(8,3) = EDXmassvar(2)./EDXmassvar(3);
            %                         ERRcoeffmatrix(8,3) = sqrt(((EDXmassvar(2).*0.05)./EDXmassvar(3)).^2 + (EDXmassvar(2).*EDXmassvar(3).*0.05./EDXmassvar(3).^2).^2);
            %                         ansmatrix(8) = fsolutions(2);
            %                     else %No O
            %                         %systemout(8) = f(3);
            %                         coeffmatrix(8,3) = 1;
            %                         ERRcoeffmatrix(8,3) = 0;
            %                         ansmatrix(8) = 0;
            %                     end
            %                 else %No N
            %                     %systemout(7) = f(2);
            %                     coeffmatrix(7,2) = 1;
            %                     ERRcoeffmatrix(7,2) = 0;
            %                     ansmatrix(7) = 0;
            %                     if ~isempty(eleidx2(eleidx2 == 3)); %O present
            %                         %systemout(8) = -f(1) + f(3) .*(mass(1)./mass(3));
            %                         coeffmatrix(8,3) = EDXmassvar(1)./EDXmassvar(3);
            %                         ERRcoeffmatrix(8,3) = sqrt(((EDXmassvar(1).*0.05)./EDXmassvar(3)).^2 + (EDXmassvar(1).*EDXmassvar(3).*0.05./EDXmassvar(3).^2).^2);
            %                         ansmatrix(8) = fsolutions(1);
            %                     else %No O
            %                         %systemout(8) = f(3);
            %                         coeffmatrix(8,3) = 1;
            %                         ERRcoeffmatrix(8,3) = 0;
            %                         ansmatrix(8) = 0;
            %                     end
            %                 end
            %             else %No C
            %                 %systemout(7) = f(1);
            %                 coeffmatrix(7,1) = 1;
            %                 ERRcoeffmatrix(7,1) = 0;
            %                 ansmatrix(7) = 0;
            %                 if ~isempty(eleidx2(eleidx2 == 2)) && ~isempty(eleidx2(eleidx2==3)); %N AND O present
            %                     %systemout(8) = -f(2) + f(3) .*(mass(2)./mass(3));
            %                     coeffmatrix(8,3) = EDXmassvar(2)./EDXmassvar(3);
            %                     ERRcoeffmatrix(8,3) = sqrt(((EDXmassvar(2).*0.05)./EDXmassvar(3)).^2 + (EDXmassvar(2).*EDXmassvar(3).*0.05./EDXmassvar(3).^2).^2);
            %                     ansmatrix(8) = fsolutions(2);
            %                 else %Either N or O or neither is present.  This doesn't check the case where only N or only O is present, but then a system with 2 mass parameters cannot be made.
            %                     %systemout(8) = f(2);
            %                     coeffmatrix(8,2) = 1;
            %                     ERRcoeffmatrix(8,2) = 0;
            %                     ansmatrix(8) = 0;
            %                 end
            %             end
            %
            %             semeleidx = eleidx2(eleidx2>3); %this gets the indexes of all non-zero elements that were analyzed with SEM
            %             if isempty(semeleidx)
            %                 semeleidx = 0;
            %             end
            %
            %             for yy = 1:(length(semeleidx)-1);
            %                 %systemout(8+i) = -f(semeleidx(i)) + f(semeleidx(i+1)).*(mass(semeleidx(i))./mass(semeleidx(i+1)));
            %                 if 8+yy == 100 %14 %14 when using the critically determined system
            %                     break
            %                 else
            %                 coeffmatrix(8+yy,semeleidx(yy+1)) = EDXmassvar(semeleidx(yy))./EDXmassvar(semeleidx(yy+1));
            %                 ERRcoeffmatrix(8+yy,semeleidx(yy+1)) = sqrt(((EDXmassvar(semeleidx(yy)).*0.05)./EDXmassvar(semeleidx(yy+1))).^2 + (EDXmassvar(semeleidx(yy)).*EDXmassvar(semeleidx(yy+1)).*0.05./EDXmassvar(semeleidx(yy+1)).^2).^2);
            %                 ansmatrix(8+yy) = fsolutions(semeleidx(yy));
            %                 end
            %             end
            %
            %             numeqns_sem = 8+(length(semeleidx)-1);
            %             sem_eleidx_zeros = eleidx_zeros(eleidx_zeros > 3);
            %             for yy = 1:length(sem_eleidx_zeros) %this explicitly makes f of any zero element, equal to zero
            %                 %systemout(numeqns_sem+j) = f(eleidx_zeros(j));
            %                 if  numeqns_sem+yy == 100 %14 %14 when using the critically determined system
            %                     break
            %                 else
            %                 coeffmatrix(numeqns_sem+yy,sem_eleidx_zeros(yy)) = 1;
            %                 ERRcoeffmatrix(numeqns_sem+yy,sem_eleidx_zeros(yy)) = 0;
            %                 ansmatrix(numeqns_sem+yy) = 0;
            %                 end
            %             end
            %             tempERRansmatrix = zeros(size(ansmatrix));
            %             tempERRansmatrix(1:6) = ERRansmatrix(1:6);
            %             ERRansmatrix = tempERRansmatrix;
            %             ERRansmatrix(7:end) = ansmatrix(7:end).*0.05;
            %             %numeqns_sem_zeros = numeqns_sem+length(eleidx_zeros);
            %             %systemout(numeqns_sem_zeros+1) = 1- sum(f(1:14));
            %             finaleqpos = size(coeffmatrix,1);
            %             coeffmatrix(finaleqpos+1,:) = ones(1,14);
            %             ERRcoeffmatrix(finaleqpos+1,:) = 0;
            %             ansmatrix(finaleqpos+1) = 1;
            %             ERRansmatrix(finaleqpos+1) = 0;
            % [coeffmatrix] * x = [ansmatrix]
            % x = [coeffmatrix]^-1 * [ansmatrix]
            
            %you cannot invert a non-square matrix.  In this code usually
            %more equations are produced than unknowns.  14 unknowns are
            %considered and so the coefficient matrix must be pared down
            %accordingly
            
            %             for oo = 1:size(coeffmatrix,1);
            %                 if size(coeffmatrix,1) == 14
            %                     break
            %                 elseif sum(coeffmatrix(oo,:)) == 1
            %                     coeffmatrix = cat(1,coeffmatrix(1:oo-1,:),coeffmatrix(oo+1:end,:));
            %                     ERRcoeffmatrix = cat(1,ERRcoeffmatrix(1:oo-1,:),ERRcoeffmatrix(oo+1:end,:));
            %                     ansmatrix = cat(1,ansmatrix(1:oo-1,:),ansmatrix(oo+1:end,:));
            %                     ERRansmatrix = cat(1,ERRansmatrix(1:oo-1,:),ERRansmatrix(oo+1:end,:));
            %                     oo = oo - 1;
            %                 end
            %             end
            
            %Since coeffmatrix' * coeffmatrix is invertible, the
            %moore-penrose pseudoinverse (used to solve an overdetermined
            %matrix equation) is calculated by (A' * A)^-1 * A'
            %MPinvcoeffmatrix =pinv(coeffmatrix);
            
            %invcoeffmatrix = inv(coeffmatrix);
            %fmatrix = ansmatrix\coeffmatrix; %fmatrix should be the exact same as fsolutions
            
            %from
            %www.cs.cornell.edu/courses/cs3220/2008su/slides/linearsystemerror.pdf
            %relerrorfsolutions = norm(coeffmatrix).*norm(invcoeffmatrix).*((norm(ERRansmatrix)./norm(ansmatrix)) + (norm(ERRcoeffmatrix)./norm(coeffmatrix)));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%End of Solver Stuff%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            partname = ['particle',num2str(j)];
                        currelefracmap = zeros(EDXrow,EDXcol);
            for g = 1:14
                currpartfraclindex = find(SEMlabelmap==j);
                Elefracmap(currpartfraclindex + (g-1).*EDXrow.*EDXcol) = fsolutions(g);    
            end
            maskedRGBmat = tempRGBmat.*partmask;
            sp2idx = find(maskedRGBmat==1);
            orgidx = find(maskedRGBmat==2);
            inorgidx = find(maskedRGBmat==3);
            
            
            %cropping particles from stxm RGB image
            croppedparts.(partname) = maskedRGBmat(cropvals(1):cropvals(2),cropvals(3):cropvals(4));
            
            %transforming and reshaping clipped SEM image into STXM coordinate space
            stxmresize_maskedRGBmat = imresize(maskedRGBmat,stxmresizeratio); %this makes a stxm image, similarly sized to an SEM image for purposes of making a reference image
            stxmresize_reference = imref2d(size(stxmresize_maskedRGBmat)); %makes an image reference to which imwarp will follow
            semclip_inv = imwarp(SEMclip,tforminv,'OutputView',stxmresize_reference);  %this applies any stretching or rotation that was identified when matching the two types of images before
            semclip_inv_revresize = imresize(semclip_inv,stxmresizeratio_reverse); %this is the final step to change an image in the SEM coordinate system into the STXM coordinate system
            croppedparts_sem.(partname) = semclip_inv_revresize(cropvals(1):cropvals(2),cropvals(3):cropvals(4)); %since we transformed SEMclip into STXM coords, we can use the cropping values found from the STXM image
            
            
            
            %sp2frac is the fraction of carbon that is attributed to soot
            if (length(sp2idx)+length(orgidx)) > 0
                sp2frac(j) = (length(sp2idx).*2)./(length(sp2idx).*2+length(orgidx).*1.3); %assuming density of soot to be 2 and density of organics to be 1.3,  this also assumes uniform thickness which is a big assumption
                
            else %this else catches the rare instance of a completely inorganic particle
                sp2frac(j) = 0;
            end
            
            %inorgfrac is the mass fraction of inorganics with respect to
            %the entire particle.
            inorgfrac(j) = (length(inorgidx).*1.77)./(length(sp2idx).*2+length(orgidx).*1.3+length(inorgidx).*1.77); %1.77 is the density of ammonium sulfate, O'Brien 2015 assumes this density for inorganics.
            
        end
        
        
        
        Crelpercentmap = Elefracmap(:,:,1);
        Convratiomap = Cmassmap./Crelpercentmap;
        nanmap = isnan(Convratiomap);
        Convratiomap(nanmap==1) = 0;
        Convratiomap(isinf(Convratiomap)==1) = 0;
        
        Elemassmap = zeros(EDXrow,EDXcol,EDXele);
        for k = 1:size(Elefracmap,3);
            Elemassmap(:,:,k) = Elefracmap(:,:,k).*Convratiomap;
        end
        
        
        %making mass.(element) struct
        numelements = length(mfracelelist);
        for qq = 1:numelements;
            evalexpression = sprintf('%s','mass.',mfracelelist{qq},'=Elemassmap(:,:,qq);');
            eval(evalexpression);
        end
        
        %calculating total mass of each pixel, for loop is needed to allow
        %for an unknown number of elements
        massfieldnames = fieldnames(mass);
        mass.tot = 0;
        for qqq = 1:numelements;
            mass.tot = mass.tot + mass.(massfieldnames{qqq});
        end
        
        
        for m = 1:numelements;
            totmass.(massfieldnames{m}) = sum(sum(mass.(massfieldnames{m}))); %u^a
        end
        totmass.tot = sum(sum(mass.tot)); %u
        
        %particle masks
        partmask = zeros(EDXrow,EDXcol,numparticles);
        parttotmass = zeros(numparticles,1);
        area = zeros(numparticles,1);
        for mm = 1:numparticles
            currpartmask = SEMlabelmap == mm;
            partmask(:,:,mm) = currpartmask;
            parttotmass(mm) = sum(sum(partmask(:,:,mm) .* mass.tot)); % u_i
            area(mm) = sum(sum(currpartmask)) .* stxmdata.(stxmfieldnames{1}).Snew.position.xstep .* stxmdata.(stxmfieldnames{1}).Snew.position.ystep; %partmask is binary, so sum(sum( will count pixels, then multiply by pixel area (in um^2)
            
        end
        AED = sqrt(4.*area./pi()); % this answer will be in um
        
        %calculating mass fractions of each element map
        for qqqq = 1:numelements
            totmfrac.(mfracelelist{qqqq}) = totmass.(massfieldnames{qqqq}) ./ totmass.tot; %p^a
        end
        
        
        partpixelnum = zeros(max(max(SEMlabelmap)),1);
        perpixel_elefrac = partpixelnum;
        temp_perpixel_elefracmap = zeros(size(SEMlabelmap));
        perpixel_elefracmap = zeros(size(Elefracmap));
        
        for l = 1:numelements;
            partelemass = zeros(numparticles,1);
            for ll = 1:numparticles
                partelemass(ll) = sum(sum(partmask(:,:,ll).*mass.(massfieldnames{l}))); %u^a_i
                
                currpartmask = SEMlabelmap == ll;
                partpixelnum(ll) = numel(find(SEMlabelmap==ll));
                perpixel_elefrac(ll) = sum(sum(Elefracmap(:,:,l).*currpartmask))./partpixelnum(ll)./partpixelnum(ll); %adds up masked elemental fraction map for a given particle, divides by number of pixels once to account for the summation and a second time to get per-pixel mfracs
                temp_perpixel_elefracmap(SEMlabelmap == ll) = perpixel_elefrac(ll);
                perpixel_elefracmap(:,:,l) = temp_perpixel_elefracmap + perpixel_elefracmap(:,:,l);
                
            end
            partmass.(massfieldnames{l}) = partelemass; %u^a_i
        end
        partmass.tot = parttotmass; % u_i
        
        partmfrac = zeros(numparticles,1);
        for w = 1:numparticles
            partmfrac(w) = partmass.tot(w) ./ totmass.tot; %p_i
        end
        
        partelemfrac = zeros(numparticles,1);
        elemfracmatrix = zeros(numparticles,14);
        for p = 1:numelements;
            for pp = 1:numparticles;
                partelemfrac(pp) = partmass.(massfieldnames{p})(pp) ./ partmass.tot(pp); %p^a_i
                elemfracmatrix(:,p) = partelemfrac;
            end
            mfrac.(massfieldnames{p}) = partelemfrac; %p^a_i
            
        end
        
        
        
        %Calculating entropies
        Hi = zeros(numparticles,1);
        Hy = 0;
        for ww = 1:numelements;
            tempHi = (-mfrac.(massfieldnames{ww}) .* log(mfrac.(massfieldnames{ww})));
            tempHi(isnan(tempHi) == 1) = 0;
            Hi = Hi + tempHi;
            Di = exp(Hi);
            
            tempHy = (-totmfrac.(massfieldnames{ww}) .* log(totmfrac.(massfieldnames{ww}))); %summing over elements
            tempHy(isnan(tempHy) == 1) = 0;
            Hy = Hy + tempHy;
        end
        
        
        Ha = sum(partmfrac .* Hi);
        Da = exp(Ha);
        Dy = exp(Hy);
        
        
        Chi = (Da-1)./(Dy-1);
        
        
        for y = 1:length(mfracelelist);
            cmptest = strcmp(mfracelelist{y},massfieldnames);
            if ~any(cmptest) == 1
                totmfrac.(mfracelelist{y}) = 0; %this ensures that all elements looked at will have a value, even if I didn't select them during analysis
            end
        end
        
        %converting C, N, and O maps from OD to mfrac
        ODmapCNO = STXM_EDXmap(:,:,1:3);
        
        CODmap = STXM_EDXmap(:,:,1);
        NODmap = STXM_EDXmap(:,:,2);
        OODmap = STXM_EDXmap(:,:,3);

        
        for x = 1:3
            ODsum = sum(sum(STXM_EDXmap(:,:,x)));
            STXM_EDXmap(:,:,x) = 100.*(STXM_EDXmap(:,:,x)./ODsum).*totmfrac.C;
        end
%         for y = 1:numel(CODmap)
%                 ODsum = sum(sum(CODmap));
%                 STXM_EDXmap(:,:,1) = (CODmap./ODsum).*totmfrac.C;
%         end
%         
%         for y = 1:numel(NODmap)
%                 ODsum = sum(sum(NODmap));
%                 STXM_EDXmap(:,:,2) = (NODmap./ODsum).*totmfrac.C;
%         end
%         
%         for y = 1:numel(OODmap)
%                 ODsum = sum(sum(OODmap));
%                 STXM_EDXmap(:,:,3) = (OODmap./ODsum).*totmfrac.C;
%         end
        
        Dimap = zeros(size(SEMlabelmap));
        for q = 1:numparticles
            Dimap(SEMlabelmap==q) = Di(q);
        end
        
        %         for bb = 1:size(elemfracmatrix,1);
        %             clnaratio = elemfracmatrix(bb,8)./elemfracmatrix(bb,4);
        %             if isnan(clnaratio)
        %                 elemfracmatrix(bb,15) = 0;
        %             elseif isinf(clnaratio)
        %                 elemfracmatrix(bb,15) = 100;
        %             else
        %                 elemfracmatrix(bb,15) = clnaratio; %to probe Cl depletion, the Cl/Na ratio is calculated.  Na is on bottom because it doesn't change whereas Cl fraction will change with chemistry
        %             end
        %         end
        
        for bb = 1:size(elemfracmatrix,1);
            elemfracmatrix(bb,15) = AED(bb); %to probe Cl depletion, the Cl/Na ratio is calculated.  Na is on bottom because it doesn't change whereas Cl fraction will change with chemistry
            elemfracmatrix(bb,16) = Di(bb);
            elemfracmatrix(bb,17) = sp2frac(bb);
            elemfracmatrix(bb,18) = inorgfrac(bb);
            elemfracmatrix(bb,19) = partmass.tot(bb); %this adds the absolute mass of the particle.  from this all masses and mass fractions can be determined (with the elemental mfrac matrix as well)
        end

        mfracmatrix_columnheader = {'C_frac','N_frac','O_frac','Na_frac','Mg_frac','P_frac','S_frac','Cl_frac','K_frac','Ca_frac','Mn_frac','Fe_frac','Ni_frac','Zn_frac','AED','Di','sp2fraction','Inorganic Fraction','Particle mfrac(Pi)'};
                
        Mixing.name = foldername;
        Mixing.Hi = Hi;
        Mixing.Ha = Ha;
        Mixing.Hy = Hy;
        Mixing.Di = Di;
        Mixing.Da = Da;
        Mixing.Dy = Dy;
        Mixing.Chi = Chi;
        Mixing.Dimap = Dimap;
        Mixing.SEMlabelmap = SEMlabelmap;
        Mixing.totmfrac = totmfrac;
        Mixing.mfrac = mfrac;
        Mixing.mfracmatrix = elemfracmatrix;
        Mixing.mfracmatrix_colheader = mfracmatrix_columnheader;
        Mixing.numparticles = numparticles;
        Mixing.AED = AED;
        Mixing.totmass = totmass;
        Mixing.ODmapCNO = ODmapCNO;
        Mixing.STXMEDXmap = STXM_EDXmap;
        Mixing.Elemassmap = Elemassmap;
        Mixing.Elefracmap = Elefracmap;
        Mixing.elelist = mfracelelist;
        Mixing.STXMCmaps = stxmdata.(stxmfieldnames{1}).Snew.RGBCompMap;
        Mixing.RawCspecavg = mean(stxmdata.(stxmfieldnames{1}).Snew.spectr(:,:,1:4),3);
        Mixing.croppedparts = croppedparts;
        Mixing.croppedparts_sem = croppedparts_sem;
        
        MixingOverview(jj,1).Mixing = Mixing;
        
        tempnumparts(jj) = numparticles;
        tempDa(jj) = Mixing.Da;
        tempDy(jj) = Mixing.Dy;
        tempCHI(jj) = Mixing.Chi;
        
        totalmfracmatrix = cat(1,totalmfracmatrix,elemfracmatrix);
        totalDimatrix = cat(1,totalDimatrix,Di);
        
        clear Cmas* Con* Crel* currp* nan* Nmass* oxy* 
        clear i ii j massfield* numele* q qq qqq ss w ww x elem* mass mm
        
        
        waitbar(jj/length(filedirs));
    end
    close(hwait);
    
    numparts = sum(tempnumparts);
	DaStats = SimpleStats(tempDa);
	DyStats = SimpleStats(tempDy);
	ChiStats = SimpleStats(tempCHI);
    MixStateStats = struct('numparts',numparts,'DaStats',DaStats,'DyStats',DyStats,'ChiStats',ChiStats);

    
    for vv = 1:length(mfracelelist); %preallocating and defining before it's assigned/called
        sumMfrac.(mfracelelist{vv}) = 0;
    end
    
    for v = 1:length(MixingOverview);
        for vv = 1:length(mfracelelist);
         sumMfrac.(mfracelelist{vv}) = sumMfrac.(mfracelelist{vv}) + MixingOverview(v,1).Mixing.totmfrac.(mfracelelist{vv});
        end
        
    end
    
    currtitle = datasettitle{ee};
    currtitle = currtitle{1};
    eval([currtitle '.MixingOverview = MixingOverview;']);
    eval([currtitle '.MixStateStats = MixStateStats;']);
    eval([currtitle '.sumMfrac = sumMfrac;']);
    eval([currtitle '.totalmfracmatrix = totalmfracmatrix;']);
    eval([currtitle '.numparticles = tempnumparts;']);
    eval([currtitle '.totalDimatrix = totalDimatrix;']);
end

% combomfracmatrix = [];
% for r = 1:length(datasettitle)
% combomfracmatrix = cat(1,combomfracmatrix,wrapper.(datasettitle{r}{1}).totalmfracmatrix);
% end


clear Chi* currtitle DaStats DyStats hwait mass mfrac Mix* sumM* totm*
clear v vv y temp* num* mfrace* matf* jj k kk l ll m inputa*
clear cmpt* Da Di Dy ee errm* evalexp* foldern* Ha Hi Hy dataset*
clear filedirs* curfield varfilename varfilenamecheck manualrect


UnifyingChi; %combining mixing state over different field of views to obtain a single value
saveaswrapper;
Combined = GoAmazonCombineSites(wrapper,'all');
clear wrapper
saveaswrapper;













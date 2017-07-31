function S = LoadStackRawMulti(filedir)
%function S=LoadStackRaw(filedir)
%
%Imports STXM raw data from input directoy filedir
%filedir needs to contain the STXM header file (.hdr) and the STXM data files (.xim)
%R.C. Moffet, T.R. Henn February 2009
%
%Modified by Matthew Fraund November 2015
%
%Inputs
%------
%filedir        path to STXM raw data directory
%
%Outputs
%-------
%S              structure array containing imported STXM data
%S.spectr       STXM absorption images
%S.eVenergy     Photon energies used to record images
%S.Xvalue       length of horizontal STXM image axis in m (window range)
%S.Yvalue       length of vertical STXM image axis in m (window range)
%S.position     structure of positioning values
%   {xvalues,yvalues,xcenter,ycenter,xstep,ystep,xpts,ypts}
%   all values in um except x/ypts wich are in pixels
%

cd(filedir) 

FileList=ls; %compiles all directory contents into string list which is as long as the longest entry (reason for strtrim)
lsize = size(FileList,1); %length doesn't work here because the file name may be longer than the number of files

for j = 1:lsize
    currentfile = strtrim(FileList(j,:)); %paring off spaces before and after file name
    hdridx = strfind(currentfile,'hdr'); %finding hdr files
    if any(hdridx)
        [S.eVenergy,S.Xvalue,S.Yvalue,multiregion,S.position]=ReadHdrMulti(currentfile); %running modified ReadHdr
        break
    end
end

if multiregion == 0; %normal image
    spccnt=1;
    for i=1:lsize
        currentfile = strtrim(FileList(i,:));
        stridx=strfind(currentfile,'xim');
        if ~isempty(stridx)
            S.spectr(:,:,spccnt)=flipud(load(currentfile)); %numbers increase going downards in matricies but upwards in the XY translator, flipud fixes this
            spccnt=spccnt+1;
        end
    end
elseif multiregion == 1 %display error message
    errormsg = sprintf('%s',filedir,' is a multistack dir, run stxmsort with multihdrsplit first');
    errordlg(errormsg);
end
            
        
% truncate crashed stacks:
if size(S.spectr,3)<length(S.eVenergy)
    S.eVenergy((size(S.spectr,3)+1):length(S.eVenergy))=[];
end

end
function [S,Snew,Mixing,Particles] = LoadImageRawMixingStateOutputCNO(filedir,names)
%function [S,Snew,Mixing,Particles] = LoadImageRawGuiMixingStateOutput(filedir,names)
%
%Imports STXM raw data from input directoy filedir
%filedir needs to contain the STXM header file (.hdr) and the STXM data files (.xim)
%R.C. Moffet, T.R. Henn February 2009
%slightly modified by MWF 2015
%
%Inputs
%------
%filedir        path to STXM raw data directory
%names           list of filenames as a 1xn or nx1 cell array of strings
%
%Outputs
%-------
%S              structure array containing imported STXM data
%S.spectr       STXM absorption images
%S.eVenergy     Photon energies used to record images
%S.Xvalue       length of horizontal STXM image axis in m
%S.Yvalue       length of vertical STXM image axis in m

cd(filedir)

%FileStruct=dir;
%numobj=length(names);
cnt=1;
for i = 1:length(names) %% loops through stack folders in raw data folder
    bidx=strfind(names{i},'.hdr');
    if ~isempty(bidx)
        NameBase{cnt}=names{i}(1:end-4); %loop is small (~4-8 iterations), preallocating takes longer
        cnt=cnt+1;
    end
end
S.eVenergy=zeros(1,length(NameBase));
for i=1:length(NameBase)
    if i==1
        temp=flipud(load(sprintf('%s_a.xim',NameBase{i})));
        S.spectr=zeros(length(temp(:,1)),length(temp(1,:)),length(NameBase));
        S.spectr(:,:,i)=temp;
        [S.eVenergy(i),S.Xvalue,S.Yvalue,~,S.position]=ReadHdrMulti(sprintf('%s.hdr',NameBase{i}));
        S.particle=NameBase{i};
    else
        S.spectr(:,:,i)=flipud(load(sprintf('%s_a.xim',NameBase{i})));
        [S.eVenergy(i),S.Xvalue,S.Yvalue,~,S.position]=ReadHdrMulti(sprintf('%s.hdr',NameBase{i}));
    end
end
[S.eVenergy,idx]=sort(S.eVenergy);
S.spectr=S.spectr(:,:,idx);
S=AlignStack(S);
if length(S.eVenergy)<5
    Snew=OdStack(S,'map',0);
else
    Snew=OdStack(S,'O',0);
end
%             cd(FinDir)
% load(sprintf('%s%s','F',S.particle));
Snew=CarbonMapsSuppFigs(Snew);
Snew = CNOeleMaps(Snew);
%STACKLab(Snew)% spccnt=1;
[Mixing,Particles] = MixingStateCNO(Snew);
try
    S.binmap=Snew.binmap;
catch
    Snew=Snew;
end
save(sprintf('%s%s','F',S.particle))

% for j=1:length(NameBase)
%     S.particle=NameBase{j};
%     for i=3:length(FileStruct)
%         stridx=findstr(FileStruct(i).name,sprintf('%s_a.xim',NameBase{j}));
%         hdridx=findstr(FileStruct(i).name,sprintf('%s.hdr',NameBase{j}));
%
%         if ~isempty(stridx)
%             S.spectr(:,:,j)=flipud(load(FileStruct(i).name));
%             spccnt=spccnt+1;
%         elseif ~isempty(hdridx)
%             [S.eVenergy(j),S.Xvalue,S.Yvalue]=ReadHdr(FileStruct(i).name);
%         end
%     end
% end
% xAxislabel=[0,S.Xvalue];
% yAxislabel=[0,S.Yvalue];
% figure,
% imagesc(xAxislabel,yAxislabel,S.spectr)
% axis image
% colorbar
% title(sprintf('%s, Raw Transmission Image, %g eV',name,S.eVenergy),'Interpreter', 'none','FontSize',14,'FontWeight','normal')
% colormap gray
% xlabel('X-Position (m)','FontSize',14,'FontWeight','normal')
% ylabel('Y-Position (m)','FontSize',14,'FontWeight','normal')
% if figsav==1
%     filename=sprintf('%s\%s',varargin{1},name);
%     saveas(gcf,filename,'png');
% end

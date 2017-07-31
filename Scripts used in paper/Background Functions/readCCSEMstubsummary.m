function [pixelsize,startdate] = readCCSEMstubsummary(txtfile)

%this reads the stubsummary and finds the pixelsize in um

fid = fopen(txtfile,'r');
rowcnt = 1;

while feof(fid) == 0;
    line = fgets(fid); %fgets advances the line of text being read
    
    pixidx = strfind(line,'Pixel size');
    starttimeidx = strfind(line,'Starting Time');
    
    if ~isempty(pixidx)
        strpixsize = line(pixidx+24:pixidx+27);
        pixelsize = str2double(strpixsize);
    elseif ~isempty(starttimeidx)
        startdate_str = line(starttimeidx+28,starttimeidx+37);
        startdate = datetime(startdate_str,'InputFormat','MM-dd-yyyy');
    end
end
fclose(fid);


end
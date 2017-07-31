function [data,colheaders,textdata] = readCCSEMcsv(csvfile)

fid = fopen(csvfile,'r');
rowcnt = 1;
hdrflag = 0;
datacnt = 1;
textdata = cell(14,1); %textdata will be at least 14 rows, more if additional phases were used
while feof(fid) == 0;
    line = fgets(fid); %fgets advances the line of text being read
    
    if hdrflag == 1;
        rowcnt = rowcnt + 1;
        continue
    end
    
    headertest = strfind(line,'Part#');
    if ~isempty(headertest)
        hdrrow = rowcnt;
        colheaders = strsplit(line,',');
        hdrflag = 1;
        continue
    end
    
    textdata{rowcnt,1} = line;
    
    rowcnt = rowcnt + 1;
end

%data = cell2mat(celldata);
fclose(fid);

colnum = length(colheaders)-1;

rowoffset1 = hdrrow;
rownum = rowcnt - 1;
data = csvread(csvfile,rowoffset1,0,[rowoffset1 0 rownum colnum]);



end
function [fileExist, probeID, data, nSamples] = gseDataRead(MatrixFileFullName);

fid = fopen(MatrixFileFullName);
fid
if (fid < 0)
    data = [];
    fileExist = 0;
    probeID = []; 
    nSamples = 0;
else
    fileExist = 1;
    tmpL = [];
    while length(strmatch('!series_matrix_table_begin', tmpL)) == 0
        tmpL = fgetl(fid);
    end;
    tmpL = fgetl(fid); A = textscan(tmpL, '%s', 'delimiter', '\t'); 
    A = A{1};
    nSamples = length(A) - 1;
    textStr = strcat('%s ', repmat('%f ', 1, nSamples));
    A = textscan(fid, textStr, 'delimiter', '\t');
    fclose(fid);
    
    probeID = A{1};
    nSamples
    length(probeID)
    data = [];
    for i = 1 : nSamples
        data = [data A{i+1}];
    end;
    
    if (length(strmatch('!series_matrix_table_end', probeID{end})) > 0)
        probeID = probeID(1 : end-1);
        data = data(1:end-1, :);
    end;
    
end;
    
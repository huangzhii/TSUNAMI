function [fileExist, probeID, data, nSamples, properties] = gseFileRead(MatrixFileFullName, PropertyList);
        properties = {};    

fid = fopen(MatrixFileFullName);
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

        if (length(PropertyList) > 0)
            for i = 1 : length(PropertyList)
                if (length(strmatch(PropertyList{i}, tmpL)) > 0)
                    A = textscan(tmpL, '%s', 'delimiter', '\t', 'bufSize', 4096*4-1); A = A{1};
                    properties = [properties, {A}]; 
                end;
            end;
        end;
    end;
    tmpL = fgetl(fid); A = textscan(tmpL, '%s', 'delimiter', '\t'); 
    A = A{1};
    nSamples = length(A) - 1;
    textStr = strcat('%s ', repmat('%f ', 1, nSamples));
    A = textscan(fid, textStr, 'delimiter', '\t');
    fclose(fid);
    
    probeID = A{1};
    for i = 1 : length(probeID)
        probeID{i} = probeID{i}(2:end-1);
    end;
    data = [];
    for i = 1 : nSamples
        data = [data A{i+1}];
    end;
    
    if (length(strmatch('!series_matrix_table_end', probeID{end})) > 0)
        probeID = probeID(1 : end-1);
        data = data(1:end-1, :);
    end;
    
end;
    
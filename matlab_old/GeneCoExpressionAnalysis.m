% This code will scan the sorted unique RefGene table for a fixed group of genes (bin) in each
% chrm in slding window, then check their correlation in a gene exp microarray. Output the
% first gene name, chr, and the ratio of pairs (within that group) with PCC > cutoff
% The NaN entries will be removed in the last step, and PCC of a gene pair with >5
% sample values will be computed.
% Change file name, path, GPL, bin, cutoff and filtering percentile,output file name.
% First version in May 1st,2013.

close all;
clear all;

% T = readtable('fpkm_table.csv','ReadVariableNames',1,'ReadRowNames',1);
% geneExpression = table2array(T);
% 
% T = readtable('rows-genes.csv','ReadVariableNames',1);
% geneSym = table2cell(T(:,4));

load('rnaData.mat');
geoFile = 'Stomach_Cancer.txt';


%Filter out probes w/o gene names, low exp value and low variance data.
%geoDataSort('---',:)= [];%remove data with no gene names
%remove data with lowest 20% absolute exp value shared by all samples
[mask, geoDataFilter, geneSymFilter]= genelowvalfilter(transpose(rna),geneId,'percentile',20);
%remove data with lowest 10% variance across samples
[mask2, geoDataFilter2, geneSymFilter2] = genevarfilter(geoDataFilter,geneSymFilter);

expData = double(geoDataFilter2);


%remove multiple gene exp values, retain only the gene exp value with
%highest mean for each gene using code HighExpressionProbes.m
[ind1, uniGene] = HighExpressionProbes(geneSymFilter2, geneSymFilter2, expData);
tmpExp = expData(ind1,:); %rows re-sorted to the alphabetic order of uniGene.


nSample = size(tmpExp, 2);

[sortMean, sortInd] = sort(mean(tmpExp, 2), 'descend');

topN = min(2000,size(tmpExp,1));

finalExp = tmpExp(sortInd(1:topN), :); % sorted expression
finalSym = uniGene(sortInd(1:topN));

tic

cMatrix = massivePCC_withoutNaN(finalExp); % calculating Pearson Correlation Coeff
%cMatrix = corr(transpose(finalExp),'type','Spearman');
cMatrix(1 : size(cMatrix,1)+1 : end) = 0; % make the diagonal be 0

toc
% tic
% Y = mdscale(1-abs(cMatrix), 3);
% toc

%%%%%%% if weigth normalization is needed, use below %%%%%%%%%%%%%%%%%%%%
% cMatrix = abs(cMatrix);
% D = sum(cMatrix);
% D_half = 1./sqrt (D);
% for i = 1 : size(cMatrix, 1)
%     cMatrix(i, :) = cMatrix(i, :) * D_half(i);
% end;
% for i = 1 : size(cMatrix, 1)
%     cMatrix(:, i) = cMatrix(:, i) * D_half(i);
% end;

%%

%%%%%%%% Step 2 - identify co-expression modules %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Algorithm parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for gamma = 0.5:0.05:0.5
% for gamma = 0.5:0.05:0.85
% for gamma = 0.2:0.05:0.45

    tic
    t = 1; lambda = 1;
    %%%%%%%% Run the algorithm
    C = localMaximumQCM(abs(cMatrix), gamma, t, lambda);
    toc
    
%     C =
% 
%   1×1768 cell array
% 
%   Columns 1 through 7
% 
%     {1×37 double}    {1×14 double}    {1×22 double}    {1×7 double}    {1×19 double}    {1×22 double}    {1×9 double}
% 
%   Columns 8 through 14
% 
%     {1×8 double}    {1×13 double}    {1×4 double}    {1×15 double}    {1×29 double}    {1×8 double}    {1×12 double}

    %%%%%%%% Step 3 - Merge the overlapped networks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%% Allowed overlap ratio threshold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta = 0.4; minClusterSize = 10;

    %%%%%%%% Sort and iteratively merge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sizeC = zeros(1, length(C)); % 1×1768
    for i = 1 : length(C)
        sizeC(i) = length(C{i});
    end;

    [sortC, sortInd] = sort(sizeC, 'descend');
    C = C(sortInd); % Still C, but sorted based on number of elements in each cell

    ind = find(sortC >= minClusterSize);

    mergedCluster = C(ind);
    mergeOccur = 1;
    currentInd = 0;

    disp('start merge')
    while mergeOccur == 1
        mergeOccur = 0;
        while currentInd < length(mergedCluster)
            currentInd = currentInd + 1;
            excludeInd = [];
            if (currentInd < length(mergedCluster))
                keepInd = 1 : currentInd;
                for j = currentInd+1 : length(mergedCluster)
                    interCluster = intersect(mergedCluster{currentInd}, mergedCluster{j});
                    if length(interCluster) >= beta*min(length(mergedCluster{j}), length(mergedCluster{currentInd}))
                        mergedCluster{currentInd} = union(mergedCluster{currentInd}, mergedCluster{j});
                        mergeOccur = 1;
                    else
                        keepInd = [keepInd, j];
                    end;
                end;
                mergedCluster = mergedCluster(keepInd);
                length(mergedCluster);
            end;
        end;
        sizeMergedCluster = zeros(1, length(mergedCluster));
        for i = 1 : length(mergedCluster)
            sizeMergedCluster(i) = length(mergedCluster{i});
        end;
        [sortSize, sortMergedInd] = sort(sizeMergedCluster, 'descend');
        mergedCluster = mergedCluster(sortMergedInd);
        currentInd = 0;
    end;
    toc

    filepath = './';
    result  = geoFile(1:length(geoFile)-4);
    outputName = strcat(filepath, result, '_CoexpressedCluster_unNormalized', num2str(round(gamma*100)),'_',num2str(minClusterSize),'.txt');
    fid = fopen(outputName, 'w');
    coexpression_name = cell(1,length(mergedCluster));
    coexpression_value = cell(1,length(mergedCluster));

    for i = 1 : length(mergedCluster)
        tmp = mergedCluster{i};
        fprintf(fid,'%s\t%d\t',strcat('cluster ',num2str(i)),length(tmp));
        for j = 1 : length(tmp)
            fprintf(fid, '%s\t', finalSym{tmp(j)});
            coexpression_name{i}{j}= finalSym{tmp(j)};
        end;
        fprintf(fid, '\r\n');
        coexpression_value{i} = finalExp(tmp,:);


    end;

    Matrix = [];
    for i = 1: length(mergedCluster)

        X = coexpression_value{i};
        mu = mean(X,2);
        stddev = std(X,0,2);
        XNorm = bsxfun(@minus,X,mu); % subtract mean of each feature from original value
        XNorm = bsxfun(@rdivide,XNorm,stddev); % divide by standard deviation

        [U,S,V] = svd(XNorm);
        Matrix = [Matrix;transpose(V(:,1))];
    end

    outputMat = strcat(filepath,result, '_', num2str(round(gamma*100)),'_',num2str(minClusterSize),'_unNormalized.mat');
    save(outputMat,'coexpression_name','coexpression_value','finalExp','finalSym','Matrix');
    fclose(fid);
    
end


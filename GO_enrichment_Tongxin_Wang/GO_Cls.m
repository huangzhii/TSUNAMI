% Gene Ontology Enrichment Analysis
% deal with a cluster file with multiplr clusters
%
% GOenrichmentSaveFile: output file name
% clsInfo: info of each gene cluster
% clsGenes: gene symbols of each gene cluster e.g 40*1 cell, each cell contains some number of gene symbols
% strGene: str.uniSym, unique gene symbols in the platform
% GO: loaded GO in matlab (GO = geneont('live',true);)
% GOtotalGenesCount: pre computed gene count for each GO item in this platform
% p_thr: p value cutoff
function numGOs = GO_Cls(GOenrichmentSaveFile, clsInfo, clsGenes, strGene,...
    maxNumGO, GO, GOrelatives, GOtotalGenesCount, p_thr)
%p_thr = 0.05;
numGOs = [];
totalGO = GO.Terms(end).id;
fid = fopen(GOenrichmentSaveFile, 'w');
% analyse each cluster
for i = 1:1:length(clsGenes)
    %fprintf('Cluster %d/%d\n', i, length(clsGenes));
    
    %%%do enrichment analysis
    [sortedpvalues, listGenesCount, totalGenesCount, goIdx] =...
        enrichment_GO(strGene.uniSym, clsGenes{i}, GOrelatives, GOtotalGenesCount);

    
    %old way
%     [sortedpvalues, listGenesCount, totalGenesCount, goIdx] =...
%         enrichment_GO_v1(strGene.uniSym, clsGenes{i}, GOrelatives, totalGO);

    
    %p-value correction here
    adjustededpvalues = mafdr(sortedpvalues,'BHFDR',true);
    
    numGO = min(maxNumGO, sum(adjustededpvalues <= p_thr));
    numGOs = [numGOs; numGO];
    
    %%% output
    %output cluster info
    clsInfotmp = clsInfo{i};
    for j = 1:1:length(clsInfo{i})
        fprintf(fid, '%s\t', clsInfotmp{j});
    end
    fprintf(fid, '\n');
    % output GO enrichment result
    for j = 1:numGO
        term = goIdx(j);
        if (isempty(GO(term).Term))
            continue;
        end
        fprintf(fid, '%s\t%-1.5d\t%-1.5d\t%d/%d\t%s\n',...
                char(num2goid(term)), sortedpvalues(j), adjustededpvalues(j),...
                listGenesCount(j), totalGenesCount(j),...
                GO(term).Term.name);
    end
end
fclose(fid);

end

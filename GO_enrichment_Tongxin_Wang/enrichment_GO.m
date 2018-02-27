%enrichment analysis of gene ontology
%geneSym is the unique gene symbols of the whole platform e.g. 10000*1 cell
%geneList, list of genes for enrichment analysis, e.g. 47*1 cell
%GOrelatives pre computed GOrelatives result
%GOtotalGenesCount pre computed gene count for each GO item in this platform
function [sortedpvalues, listGenesCount, totalGenesCount, goIdx] = ...
    enrichment_GO(geneSym, geneList, GOrelatives, GOtotalGenesCount)
totalGenesCount = GOtotalGenesCount;
totalGO = length(totalGenesCount);
listGenesCount  = zeros(totalGO,1);  % a vector of GO term counts for up-regulated genes.

for i = 1:length(geneList)
    if isKey(GOrelatives,geneList{i})
        goid = GOrelatives(geneList{i});
        listGenesCount(goid) = listGenesCount(goid) + 1;
    end
end

goIdx = find(totalGenesCount > 0);
totalGenesCount = totalGenesCount(goIdx);
listGenesCount = listGenesCount(goIdx);

% result on the matlab page seems wrong...
% gopvalues = hygepdf(listGenesCount,max(totalGenesCount),...
%                         max(listGenesCount),totalGenesCount);
gopvalues = hygepdf(listGenesCount,length(geneSym),...
                        totalGenesCount,length(geneList));



[sortedpvalues, goIdx_tmp] = sort(gopvalues);
listGenesCount = listGenesCount(goIdx_tmp);
totalGenesCount = totalGenesCount(goIdx_tmp);
goIdx = goIdx(goIdx_tmp);

% only output GO items that have been enriched
idx_non0 = find(listGenesCount > 0);
sortedpvalues = sortedpvalues(idx_non0);
listGenesCount = listGenesCount(idx_non0);
totalGenesCount = totalGenesCount(idx_non0);
goIdx = goIdx(idx_non0);

end



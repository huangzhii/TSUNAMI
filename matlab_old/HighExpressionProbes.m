function [ind, uniGene] = HighExpressionProbes(Genes, Probes, Data);
meanData = mean(Data, 2);
[sGenes, sInd] = sort(Genes);
sProbes = Probes(sInd);
sMean = meanData(sInd);
[uniGene, uniInd] = unique(sGenes);
uniInd = [0; uniInd];
tmpInd = zeros(1, length(uniGene));
for i = 1 : length(uniInd)-1
    [maxV, maxInd] = max(sMean(uniInd(i)+1:uniInd(i+1)));
    tmpInd(i) = uniInd(i) + maxInd;
end;
% [maxV, maxInd] = max(sMean(uniInd(end):end));
% tmpInd(end) = uniInd(end) + maxInd - 1;
ind = sInd(tmpInd);


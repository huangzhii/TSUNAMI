% Gene Ontology Enrichment Analysis
% deal with multiple cluster files in a folder
%
% p_thr: p value cutoff
function fun_GO(folder, geoID, subfolder, mat_file_end, GPLFile, cls_file_end, p_thr)
%%%%%%%%%%%%%%%%%%%%%
maxNumGO = 10;
%p_thr = 0.05;
%%%%%%%%%%%%%%%%%%%%%
%str.finalExp, str.finalSym, str.uniExp, str.uniSym
strGene = load(fullfile(folder, geoID,[geoID,mat_file_end]));
%str.GOrelatives strGO.totalGenesCount
strGO = load(fullfile(folder,geoID,[GPLFile(1:end-4),'_GOrelatives.mat']));
%construct gene Symbol and GOid mapping
GO = geneont('live',true);
%%%%%%%%%%%%%%%%%%%%%
% folders in subfolder
% folder contains multiple cluster files

geneClsFiles = dir([fullfile(folder, geoID, subfolder),'/*',cls_file_end]);
    
%%%%%%%%%%
for i=1:1:length(geneClsFiles)
    mf_GOSaveFile = fullfile(geneClsFiles(i).folder,[geneClsFiles(i).name(1:end-length(cls_file_end)),'.GOmf']);
    bp_GOSaveFile = fullfile(geneClsFiles(i).folder,[geneClsFiles(i).name(1:end-length(cls_file_end)),'.GObp']);
    cc_GOSaveFile = fullfile(geneClsFiles(i).folder,[geneClsFiles(i).name(1:end-length(cls_file_end)),'.GOcc']);
    
    clsFileName = fullfile(geneClsFiles(i).folder,geneClsFiles(i).name);
    %determine if the file is cluster of genes or cluster of eigenGenes
    if isEigenGeneFile(clsFileName) == 1
        [clsInfo, clsGenes] = eigenGeneClsRead(clsFileName);
        for j=1:1:length(clsGenes)
            clsGenes{j} = vertcat(clsGenes{j}{:});
        end
    elseif isEigenGeneFile(clsFileName) == 0
        [clsInfo, clsGenes] = geneClsRead(clsFileName);
    else
        printf("file read error\n")
        continue;
    end
    
    %GO Enrichment Analysis for each file
    fprintf('Doing %s\n',fullfile(geneClsFiles(i).name));
    %molecular function
    fprintf('molecular function\n');
    mf_numGOs = GO_Cls(mf_GOSaveFile, clsInfo, clsGenes, strGene, maxNumGO, GO, ...
        strGO.mf_GOrelatives, strGO.mf_totalGenesCount, p_thr);
    %biological process
    fprintf('biological process\n');
    bp_numGOs = GO_Cls(bp_GOSaveFile, clsInfo, clsGenes, strGene, maxNumGO, GO, ...
        strGO.bp_GOrelatives, strGO.bp_totalGenesCount, p_thr);
    %cellular component
    fprintf('cellular component\n');
    cc_numGOs = GO_Cls(cc_GOSaveFile, clsInfo, clsGenes, strGene, maxNumGO, GO, ...
        strGO.cc_GOrelatives, strGO.cc_totalGenesCount, p_thr);
end

end

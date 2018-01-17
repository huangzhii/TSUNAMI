% dataMat is a matrix of N-by-M containing NaN entries (N > 1)
% dataVec is a row vector of 1-by-M containing NaN entries
% PCC_vec is a column vector of N-by-1 which is the PCC between each row
% of dataMat and dataVec
% minCommonThresh is the minimum number of data pairs required for
% reporting a meaningful PCC value
function PCC_vec = massivePCC_withNaN(dataMat, dataVec, minCommonThresh);
[nRow, nCol] = size(dataMat);

non_NaN_mat_mask = 1 - isnan(dataMat);
non_NaN_vec_mask = 1 - isnan(dataVec);

common_Mask = non_NaN_mat_mask .* repmat(non_NaN_vec_mask, nRow, 1);

common_Count = sum(common_Mask, 2);

newDataMat = dataMat; % .* common_Mask; 
newDataMat(common_Mask < 0.5) = 0;
revisedMeanDataMat = sum(newDataMat, 2)./common_Count;
newDataMat = newDataMat - repmat(revisedMeanDataMat, 1, nCol);
newDataMat(common_Mask < 0.5) = 0;
rowNormMat = sqrt(sum(newDataMat .* newDataMat, 2));
newDataMat = newDataMat ./ repmat(rowNormMat, 1, nCol);

newDataVec_Mat = repmat(dataVec, nRow, 1); % .* common_Mask;
newDataVec_Mat(common_Mask < 0.5) = 0;
revisedMeanDataVec = sum(newDataVec_Mat, 2) ./ common_Count;
newDataVec_Mat = newDataVec_Mat - repmat(revisedMeanDataVec, 1, nCol);
newDataVec_Mat(common_Mask < 0.5) = 0;
rowNormVec = sqrt(sum(newDataVec_Mat .* newDataVec_Mat, 2));
newDataVec_Mat = newDataVec_Mat ./ repmat(rowNormVec, 1, nCol);

PCC_vec = sum(newDataMat .* newDataVec_Mat, 2);

ind = find(common_Count < minCommonThresh);
PCC_vec(ind) = NaN;








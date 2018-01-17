% dataMat is a matrix of N-by-M containing NaN entries (N > 1)
% dataVec is a row vector of 1-by-M containing NaN entries
% PCC_vec is a column vector of N-by-1 which is the PCC between each row
% of dataMat and dataVec
% minCommonThresh is the minimum number of data pairs required for
% reporting a meaningful PCC value (Pearson Correlation Coefficient)
function PCC_mat = massivePCC_withoutNaN(dataMat);
[nRow, nCol] = size(dataMat);

dataMat(isnan(dataMat)) = 0;


if (nRow > 1)
    
    sumX = sum(dataMat, 2); % sum each rows to a value and becomes n * 1 vector
    meanX = sumX / nCol;
    sumV = sqrt(sum(dataMat.*dataMat, 2) - (sumX .* sumX)/nCol);
    
    PCC_mat = zeros(nRow, nRow);
    
    for i = 1 : nRow
        PCC_mat(i, i+1:end) = (dataMat(i+1:end,:) * dataMat(i,:)' - sumX(i) * sumX(i+1:end) / nCol)./(sumV(i) * sumV(i+1:end));
    end;
    PCC_mat = PCC_mat + PCC_mat';
    PCC_mat(1 : nRow+1 : end) = 1;
else
    PCC_mat = corr(dataMat(1, ind), dataVec(1, ind));
end;










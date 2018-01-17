function C = localMaximumQCM(cMatrix, gamma, t, lambda);

C = {};
%%Step 1 - find the local maximal edges
[maxV, maxInd] = max(cMatrix);
maxEdges = []; maxW = [];
for i = 1 : length(maxInd)
    if maxV(i)== max(cMatrix(maxInd(i), :));
        maxEdges = [maxEdges; maxInd(i), i];
        maxW = [maxW; maxV(i)];
    end;
end;

[sortMaxV, sortMaxInd] = sort(maxW, 'descend');
sortMaxEdges = maxEdges(sortMaxInd, :);

length(sortMaxInd)

currentInit = 1; noNewInit = 0;

nodesInCluster = [];

while (currentInit <= length(sortMaxInd)) & (noNewInit == 0)
    if (sortMaxV(currentInit) < gamma * sortMaxV(1))
        noNewInit = 1;
    else
        if (ismember(sortMaxEdges(currentInit, 1), nodesInCluster)==0 & ...
                ismember(sortMaxEdges(currentInit, 2), nodesInCluster)==0)
            newCluster = sortMaxEdges(currentInit, :);
            addingMode = 1;
            currentDensity = sortMaxV(currentInit);
            nCp = 2;
            totalInd = 1 : size(cMatrix, 1);
            remainInd = setdiff(totalInd, newCluster);
            while addingMode == 1
                
                neighborWeights = sum(cMatrix(newCluster, remainInd));
                
                [maxNeighborWeight, maxNeighborInd] = max(neighborWeights);
                c_v = maxNeighborWeight/nCp;
                alphaN = 1 - 1/(2*lambda*(nCp+t));
                if (c_v >= alphaN * currentDensity)
                    newCluster = [newCluster, remainInd(maxNeighborInd)];
                    nCp = nCp+1;
                    currentDensity = (currentDensity*((nCp-1)*(nCp-2)/2)+maxNeighborWeight)/(nCp*(nCp-1)/2);
                    remainInd = setdiff(remainInd, remainInd(maxNeighborInd));
                else
                    addingMode = 0;
                end;
            end;
            nodesInCluster = [nodesInCluster, newCluster];
            C = [C, newCluster];
            
        end;
    end;
    currentInit = currentInit + 1
end;
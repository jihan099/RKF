% This function calculates the rank of the recursive Koszul flattening of
% Perm_n by using Lemma 5.5.
function accRank = getRankSymm(M, orbitmat, n)
    checkedColList = double.empty;
    strToPrint='';
    accRank = 0;
    colIndsList = orbitmat(:, 1);

    while length(checkedColList) ~= length(colIndsList)
        % Iterate checkInds
        [numToAdd, checkedColList, blockSize, M] = checkInds(M, colIndsList, checkedColList, strToPrint, orbitmat, n);
        accRank = accRank + numToAdd;
        strToPrint = strcat('Left inds to check:', string(length(setdiff(colIndsList, checkedColList))), ...
            ', accumulated rank:', string(accRank), ', last block size:', string(blockSize));
    end
end

% This function is internally used in getRankSymm.
% Given the list of indices of columns, it chooses an index that is not
% checked yet. Then it calculates the rank of the connected component
% containing the corresponding column. The calculated rank is multipled
% by n!/|H_a| (c.f. Lemma 5.5).
function [numToAdd, checkedColList, blockSize, reducedM] = checkInds(M, colIndsList, checkedColList, strToPrint, orbitmat, n)
    % Choose a column to find a block containing it
    uncheckedColList = setdiff(colIndsList, checkedColList); % The list of unchecked indices
    disp(length(colIndsList));
    disp(length(checkedColList));
    startingCol = uncheckedColList(1); % The index of the column to check

    [~, uniIndsJ, block, oriBlock] = findBlock(M, startingCol, strToPrint); % Get the connected component containing the column
    fBlock = full(block);
    rankFBlock = rank(fBlock);

    [numrowOrbitmat, ~] = size(orbitmat);
    posCheckedIndsToAdd = rem(find(ismember(orbitmat, uniIndsJ)), numrowOrbitmat);
    for i = 1:length(posCheckedIndsToAdd)
        if posCheckedIndsToAdd(i) == 0
            posCheckedIndsToAdd(i) = numrowOrbitmat;
        end
    end
    checkedColList = unique(cat(1, checkedColList, orbitmat(posCheckedIndsToAdd, 1)));
    blockSize = length(uniIndsJ);
    reducedM = M - oriBlock;

    orbit = transpose(orbitmat(colIndsList == startingCol, :));
    [~,J,~] = find(oriBlock);
    pos = find(ismember(orbit, J)); % The variable corresponds to |H_a|.
    numToAdd = rankFBlock * factorial(n) / length(pos); % Using the formula.
    return;
end

% This function is internally used in checkInds.
% Given the index of the column, it returns the connected component
% including it.
function [uniIndsI, uniIndsJ, block, oriBlock] = findBlock(M, startingCol, strToPrint)
    [I,J,V] = find(M);
    [numrowM, numcolM] = size(M);
    inds = find(J == startingCol);
    checkedI = double.empty;
    checkedJ = [startingCol];

    indsI = I(inds);
    indsJ = J(inds);
    uniIndsI = unique(indsI);
    uniIndsJ = unique(indsJ);

    while true
        if length(uniIndsI) ~= length(checkedI)
            % Pick row to check
            uncheckedI = setdiff(uniIndsI, checkedI);
            indsToAdd = find(ismember(I, uncheckedI));
            uniIndsJ = unique(cat(1, uniIndsJ, J(indsToAdd)));
            inds = cat(1, inds, indsToAdd);
            checkedI = cat(1, checkedI, uncheckedI);
        elseif length(uniIndsJ) ~= length(checkedJ)
            % Pick col to check
            uncheckedJ = setdiff(uniIndsJ, checkedJ);
            indsToAdd = find(ismember(J, uncheckedJ));
            uniIndsI = unique(cat(1, uniIndsI, I(indsToAdd)));
            inds = cat(1, inds, indsToAdd);
            checkedJ = cat(1, checkedJ, uncheckedJ);
        else
            break;
        end
        disp(strToPrint);
    end

    inds = unique(inds);
    indsI = I(inds);
    indsJ = J(inds);
    oriBlock = sparse(indsI, indsJ, V(inds), numrowM, numcolM);

    lengIndsI = length(indsI);
    newIndsI = zeros(lengIndsI, 1);
    for i = 1:lengIndsI
        newIndsI(i) = find(uniIndsI == indsI(i), 1);
    end

    lengIndsJ = length(indsJ);
    newIndsJ = zeros(lengIndsJ,1);
    for i = 1:lengIndsJ
        newIndsJ(i) = find(uniIndsJ == indsJ(i), 1);
    end

    block = sparse(newIndsI, newIndsJ, V(inds));
end
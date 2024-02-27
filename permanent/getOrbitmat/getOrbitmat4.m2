n = 4; -- The order of the permanent tensor

-- The variable indsList represents an element in the basis B_1.
-- Given the element corresponding to indsList, this function returns the index of 
-- the column in the matrix corresponds to the recursive Koszul flattening of the permanent tensor.
findColInd = (indsList,n) -> (	
	finalColInd := 0;
	currentSize := product(for i from 1 to n-2 list binomial(n,i)) * n;
	numColCurrentPsi := 0;

	for i from 1 to n-2 do (
		numColCurrentPsi = binomial(n,i);
		currentSize = sub(currentSize / numColCurrentPsi, ZZ);
		finalColInd = finalColInd + getColIndsInPsi(indsList_(i-1), n) * currentSize;
	);
	
	finalColInd = finalColInd + indsList_(n-2)_0;
	finalColInd = finalColInd + 1;
	
	return finalColInd;
)

-- Internally used in findColInd
getColIndsInPsi = (inds, n) -> (
	return position(sort subsets(n, #inds), i->(i==inds));
)

-- Given the column index of the matrix of the recursive Koszul flattening of the permanent tensor,
-- this function returns indsList that represents an element in the basis B_1.
findIndsList = (colInd,n) -> (
	colInd = colInd - 1;
	quotientInd := 0;
	indsList := {};
	currentSize := product(for i from 1 to n-2 list binomial(n,i)) * n;
	
	for i from 1 to n-2 do (
		currentSize = sub(currentSize / binomial(n,i), ZZ);
		quotientInd = colInd // currentSize;
		colInd = colInd % currentSize;
		indsList = append(indsList, (sort subsets(n, i))_quotientInd);
	);
	
	indsList = append(indsList, {colInd});
	
	return indsList;
)

-- This function performs the permutation sigma to the basis element corresponds to indsList and returns the corresponding element.
doPermu = (indsList, sigma) -> (
	outterList := {};
	innerList := 0;
	for i in indsList do(
		innerList = {};
		for j in i do (
			innerList = sort(append(innerList, sigma_j));
		);
		outterList = append(outterList, innerList);
	);
	return outterList;
)

permutationList = permutations n; -- The list of all permutations

-- Write to a file
addedInds = {};
totalNumOfInds = n * (product for i from 1 to n-2 list binomial(n,i));
totalInds = set(1..totalNumOfInds);
filename = concatenate("orbitmat", toString n, ".dat");
filename << "" << close;
openOutAppend filename;
for i from 1 to totalNumOfInds do (
	if #addedInds == totalNumOfInds then break;
	unaddedColInd = (toList(totalInds - addedInds))_0;
	indsList = findIndsList(unaddedColInd, n);
	orbit = for sigma in permutationList list findColInd(doPermu(indsList,sigma),n);
	orbit = sort(orbit);
	addedInds = join(addedInds, unique orbit);
	for j in orbit do (
		filename << j << " ";
	);
	filename << endl;
	stdio << #addedInds << " / " << totalNumOfInds << endl;
);
filename << close;

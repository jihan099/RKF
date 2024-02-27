n = 7;			 -- The order of the permanent tensor
-- As orbitmatRedundant7 is too large, we divide into 40 parts. Each part calculates 500 rows.
-- The following is the number of the part.
partnum = 13;
chunkSize = 500;

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

-- This part constructs the list of some elements in the basis. 
recursiveInds = (listOfElts, n) -> (
	outList := {};
	currentLeng := 0;
	flatIndsList := 0;
	numOfUsedInds := 0;
	numOfIndsToAdd := 0;
	lastNum := 0;
	indsToAdd := 0;
	out := 0;
	
	for indsList in listOfElts do(
		currentLeng = #indsList;
		flatIndsList = flatten indsList;
		numOfUsedInds = length unique flatIndsList;
		for numOfReusedInds from 0 to min(numOfUsedInds, currentLeng) do (
			numOfIndsToAdd = currentLeng - numOfReusedInds;
			lastNum = (numOfUsedInds - 1) + numOfIndsToAdd;
			if lastNum > n - 1 then continue;
			
			for reusedIndsList in sort subsets(numOfUsedInds, numOfReusedInds) do (
				indsToAdd = join(reusedIndsList, toList(numOfUsedInds..lastNum));
				out = insert(-2, indsToAdd, indsList);
				outList = append(outList, out);
			);
		);
	);
	
	return outList;
)
listOfElts={{{0}, {0}}, {{1}, {0}}}; -- fixed for all n
for i from 1 to n-3 do (
	listOfElts=recursiveInds(listOfElts, n);
	for i in listOfElts do print i;
);

-- Convert the elements in listOfElts to the indices of the columns of the matrix corresponding to the recursive Koszul flattening.
colIndsList = {};
for indsList in listOfElts do (
	finalColInd = findColInd(indsList, n);
	colIndsList = append(colIndsList, finalColInd);
);

permutationList = permutations n; -- The list of all permutations

-- Write to a file

filename = concatenate("orbitmatRedundant", toString n, "part", toString partnum, ".dat");
filename << "" << close;
rownumBegin = (partnum-1)*chunkSize;
rownumEnd = min(partnum*chunkSize-1, #listOfElts-1);
openOutAppend filename;
for i from rownumBegin to rownumEnd do (
	if i > #listOfElts-1 then break;
	print (i,#listOfElts-1);
	indsList = listOfElts_i;
	orbit = sort for sigma in permutationList list findColInd(doPermu(indsList,sigma),n);
	filename << i << " ";
	for j in orbit do (
		filename << j << " ";
	);
	filename << endl;
	clearOutput;
	collectGarbage();
);
filename << close;

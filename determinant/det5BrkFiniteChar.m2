-- This code computes the determinant of the recursive Koszul flattening of det5 over the field of rational numbers.

-- This function returns a list of matrices. The i-th matrix in the list corresponds to the wedge product map with respect to the i-th basis of the vector space.
-- The variable ord is the order of the tensor and dimV is the dimension of the vector space.
getPsi = (ord, dimV) -> (
	numcol := binomial(dimV, ord);
	numrow := binomial(dimV, ord+1);
	col := 0;
	outMat := 0;
	outList := {};

	for i from 0 to dimV-1 do(
		for indBeforeWedge from 0 to numcol-1 do(
			if getSign(i, indBeforeWedge, ord, dimV) == 0 then col = matrix table (numrow, 1, (r,c)->0)
			else (
				col = matrix table (numrow, 1, (r,c)->(if r==getPos(i, indBeforeWedge, ord, dimV) then getSign(i, indBeforeWedge, ord, dimV) else 0));
			);

			if indBeforeWedge == 0 then outMat = col
			else outMat = outMat | col;
		);
		outList = append(outList, outMat);
	);

	return outList;
)

-- This function is used internally in the function getPsi.
getSign = (i, posInLexList, ord, dimV) -> (
	lexList := sort subsets(dimV, ord);
	ind := lexList_posInLexList;
	for pos from 0 to ord-1 do (
		if i == ind_pos then return 0;
		if i < ind_pos then return (-1)^pos;
	);
	return (-1)^ord;
)

-- This function is used internally in the function getPsi.
getPos = (i, posInLexList, ord, dimV) -> (
	lexList := sort subsets(dimV, ord);
	ind := lexList_posInLexList;
	indAfterWedge := {};
	for pos from 0 to ord-1 do (
		if i == ind_pos then return -1;
		if i < ind_pos then (
			indAfterWedge = insert(pos, i, ind);
			break;
		);
	);
	if i > ind_(-1) then indAfterWedge = insert(ord, i, ind);

	nextLexList := sort subsets(dimV, ord+1);
	for i from 0 to #nextLexList-1 do (
		if nextLexList_i === indAfterWedge then return i;
	);
)

-- This function takes a Kronecker product of two matrices.
kronecker = (M,N) -> (
	rowM := numRows M;
	colM := numColumns M;
	rowN := numRows N;
	colN := numColumns N;
	return matrix table(rowM * rowN, colM * colN, (i,j) -> M_(i//rowN,j//colN) * N_(i%rowN, j%colN));
)

k := QQ; 													-- The base field
dimV := 5;													-- The dimension of the vector space

-- The following code is to construct the determinant tensor det5.
S := k[a_0..a_4, b_0..b_4, c_0..c_4, d_0..d_4, e_0..e_4];
M := matrix table(dimV, dimV, (i,j)->S_(dimV*i+j));
det5 := det M;												-- This corresponds to the tensor det5.

-- As we take the recursive Koszul flattning with (p_1,p_2,p_3)=(1,2,3), we take the following values.
p_1 = 1;
p_2 = 2;
p_3 = 3;

-- Construct the wedge product maps
for i from 1 to 3 do psi_i = getPsi(p_i, dimV);

-- Construct the matrix corresponding to the recursive Koszul flattening of det5
kosz = 0;
for i from 0 to dimV-1 do (
	for j from 0 to dimV-1 do (
		for l from 0 to dimV-1 do (
			if #set(i,j,l) != 3 then continue; 				-- Continue the loop if i,j,l are not distinct.
			stdio << "In the loop " << (i,j,l) << endl;		-- Print the progress
			flat_(i,j,l) = matrix table(dimV, dimV, (ii,jj)->coefficient(a_i*b_j*c_l*d_ii*e_jj, det5));
			kron_(i,j,l) = kronecker(kronecker(kronecker(psi_1_i, psi_2_j), psi_3_l), flat_(i,j,l));
			kosz = kosz + kron_(i,j,l);
		);
	);
);

-- Factorize the determinant
factor sub(det kosz, ZZ)									-- The result is 2^1600*3^25. Hence brk(det5) >= 27 when char(K) >= 5.

-- This code computes the determinant of the recursive Koszul flattening of det4-S over the field of rational numbers where S is a simple tensor.

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
dimV := 4;													-- The dimension of the vector space

-- The following code is to construct the determinant tensor det4.
S := k[x,y,z,w][a_0..a_3, b_0..b_3, c_0..c_3, d_0..d_3];
M := matrix table(dimV, dimV, (i,j)->S_(dimV*i+j));
det4 := det M;												-- This corresponds to the tensor det4.

-- For the simple tensor S = s1 (tensor) s2 (tensor) s3 (tensor) s4, we divide into four cases according to the dimension of the vector space <s1,s2,s3,s4>.
-- We express the tensor det4-S as follows in each case.
tensorCase_1 = det4 - x*a_0*b_1*c_2*d_3;					-- when dim <s1,s2,s3,s4> = 4
tensorCase_2 = det4 - a_0*b_1*c_2*(x*d_0+y*d_1+z*d_2);		-- when dim <s1,s2,s3,s4> = 3
tensorCase_3 = det4 - a_0*b_1*(x*c_0+y*c_1)*(z*d_0+w*d_1);	-- when dim <s1,s2,s3,s4> = 2
tensorCase_4 = det4 - x*a_0*b_0*c_0*d_0;					-- when dim <s1,s2,s3,s4> = 1

-- As we take the recursive Koszul flattning with (p_1,p_2)=(1,2), we take the following values.
p_1 = 1;
p_2 = 2;

-- Construct the wedge product maps
for i from 1 to 2 do psi_i = getPsi(p_i, dimV);

-- Construct the matrix corresponding to the recursive Koszul flattening of tensorCase_n.
for n from 1 to 4 do (
	koszCase_n = 0;
	for i from 0 to dimV-1 do (
		for j from 0 to dimV-1 do (
			if #set(i,j) != 2 then continue; 				-- Continue the loop if i,j are not distinct.
			flat_(i,j) = matrix table(dimV, dimV, (ii,jj)->coefficient(a_i*b_j*c_ii*d_jj, tensorCase_n));
			kron_(i,j) = kronecker(kronecker(psi_1_i, psi_2_j), flat_(i,j));
			koszCase_n = koszCase_n + kron_(i,j);
		);
	);

	stdio << "Case " << n << ": " << factor det(koszCase_n) << endl;
);
-- The result is that the determinant is a nonzero integer in the case 2,3,4.
-- In the case 1, the determinant is a scalar multiple of (x-1) * (x-2)^4 * (x-4)^4.

-- In the case 1, when x=1
stdio << "Case 1, when x=1, rank=" << rank sub(sub(koszCase_1,x=>1),k) << endl;    -- The result is 95

-- In the case 1, when x=2
stdio << "Case 1, when x=1, rank=" << rank sub(sub(koszCase_1,x=>2),k) << endl;    -- The result is 93

-- In the case 1, when x=4
stdio << "Case 1, when x=1, rank=" << rank sub(sub(koszCase_1,x=>4),k) << endl;    -- The result is 92

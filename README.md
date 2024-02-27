# Recursive Koszul flattening of determinant and permanent tensors
The purpose of the codes in this repository is to implement the arguments in the paper 'Recursive Koszul flattening of determinant and permanent tensors'. The files with the extension .m2 should be run by Macaulay2 and the files with the extension .m should be run by Matlab. Macaulay2 codes were written based on the version 1.22 and Matlab codes were written based on the version R2023a.

## Determinants
- The codes 'determinant/det4BrkFiniteChar.m2' and 'determinant/det5BrkFiniteChar.m2' are to perform the argument in Theorem 4.4.
- The code 'determinant/det4RkChar0.m2' is to perform the argument in the characteristic zero part of Theorem 4.6.
- The code 'determinant/det4RkFiniteChar.m2' is to perform the argument in the characteristic $p\geq 3$ part of Theorem 4.6.

## Permanents
- All files in the directory 'permanent' are to get the result in Theorem 5.6 using Lemma 5.5.
- To get the result for $\text{Perm}_n$, one needs to run 'permanent/RKF_Per{n}.m'.
- To run 'permanent/RKF_Per{n}.m', the files 'permanent/spkron.m', 'getRankSymm.m', and 'orbitmat{n}.dat' should be in the working directory(current folder) of Matlab.
- The files 'orbitmat{n}.dat' are contained in 'orbitmats.zip'. After unzipping it, the file 'orbitmat{n}.dat' should be moved to the working directory of Matlab, if necessary.
- Basically, 'orbitmat{n}.dat' contains a matrix whose entries are indices of columns of the matrix corresponding to the recursive Koszul flattening of $\text{Perm}_n$.
- Each row of the matrix corresponds to an orbit of an element of the basis $\mathcal{B}_1$ with respect to the $\mathfrak{S}_n$-action in the form of corresponding column indices. We allow the elements to be repeated in the orbit. Consequently, the number of columns is $n!$.

### How to get 'orbitmat{n}.dat'
- The files in 'permanent/getOrbitmat' are to make the files 'orbitmat{n}.dat'. Note that one can just use the ones in 'orbitmats.zip'.
- For $4\leq n\leq 6$, one can get 'orbitmat{n}.dat' by running 'permanent/getOrbitmat/getOrbitmat{n}.m2' which results to create the file in the current directory of Macaulay2.
- For $n=7$, to get 'orbitmat7.dat', we should get 'orbitmatRedundant7part{m}.dat' for all $1\leq m\leq 40$ first which contains some redundant rows. The first column of the matrix in 'orbitmatRedundant7part{m}.dat' represents the row indices of the whole matrix. Those are there to check if the file is created properly. Those are removed in 'orbitmat7.dat'.
- One can get 'orbitmatRedundant7part{m}.dat' by running 'permanent/getOrbitmat/getOrbitmatRedundant7part{m}.m2' which results to create the file in the current directory of Macaulay2.
- Assuming that one has 'orbitmatRedundant7part{m}.dat' for all $1\leq m\leq 40$, one can get 'orbitmat7.dat' by running 'permanent/getOrbitmat/getOrbitmat7.m' in Matlab. The files 'orbitmatRedundant7part{m}.dat' should be in the working directory of Matlab for all $1\leq m\leq 40$ when running it.
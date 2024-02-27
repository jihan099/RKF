% Recursive Koszul Flattening of Per4
n = 4;
syms a [n];
det = det(a); % Note that later, we take abs() when we construct the matrix as we consider the permanent tensor.

% Costruct psi_i.
% It corresponds to the wedge product map with respect to the i-th basis of the vector space.
N1 = nchoosek(n, 1);
N2 = nchoosek(n, 2);
N3 = nchoosek(n, 3);
psi1 = zeros(N2, N1, n);
psi2 = zeros(N3, N2, n);

% Input psi1
psi1(:,:,1) = [0,1,0,0;0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0;0,0,0,0];
psi1(:,:,2) = [-1,0,0,0;0,0,0,0;0,0,0,0;0,0,1,0;0,0,0,1;0,0,0,0];
psi1(:,:,3) = [0,0,0,0;-1,0,0,0;0,0,0,0;0,-1,0,0;0,0,0,0;0,0,0,1];
psi1(:,:,4) = [0,0,0,0;0,0,0,0;-1,0,0,0;0,0,0,0;0,-1,0,0;0,0,-1,0];

% Input psi2
psi2(:,:,1) = [0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1;0,0,0,0,0,0];
psi2(:,:,2) = [0,-1,0,0,0,0;0,0,-1,0,0,0;0,0,0,0,0,0;0,0,0,0,0,1];
psi2(:,:,3) = [1,0,0,0,0,0;0,0,0,0,0,0;0,0,-1,0,0,0;0,0,0,0,-1,0];
psi2(:,:,4) = [0,0,0,0,0,0;1,0,0,0,0,0;0,1,0,0,0,0;0,0,0,1,0,0];

% Construct the matrix correponding to the recursive Koszul flattening
N = nchoosek(4,1)*nchoosek(4,2)*n;
kosz = sparse(N,N);
flat = zeros(n,n,n,n);
for p = 1:n
    diff1 = diff(det,a(1,p));
    for q = 1:n
        if length(unique([p,q])) ~= length([p,q])
            continue
        end
        diff2 = diff(diff1,a(2,q));
        kron12 = sparse(kron(psi1(:,:,p),psi2(:,:,q)));
        for i = 1:n
            if length(unique([p,q,i])) ~= length([p,q,i])
                continue
            end
            diff3 = diff(diff2,a(3,i));
            for j = 1:n
                if length(unique([p,q,i,j])) ~= length([p,q,i,j])
                    continue
                end
                flat(i,j,p,q) = abs(diff(diff3,a(4,j)));
            end
        end
        kosz = kosz + spkron(kron12,sparse(flat(:,:,p,q)));
    end
end

% Calculate the rank
orbitmat = readmatrix('orbitmat4.dat');
tic;
rk = getRankSymm(kosz, orbitmat, 4)    % rank: 70, lower bound of brk(Per4): 8
x = toc;
x
% This function returns the Kronecker product of two sparse matrices.
function C = spkron(A, B)
   [ia, ja, va] = find(A);
   [ra, ca] = size(A);
   [ib, jb, vb] = find(B);
   [rb, cb] = size(B);
   lengA = length(ia);
   lengB = length(ib);
   ic = zeros(1, lengA*lengB);
   jc = zeros(1, lengA*lengB);
   vc = zeros(1, lengA*lengB);
   for p = 1:lengA
       for q = 1:lengB
            ic((p-1)*lengB+q) = (ia(p)-1) * rb + ib(q);
            jc((p-1)*lengB+q) = (ja(p)-1) * cb + jb(q);
            vc((p-1)*lengB+q) = va(p) * vb(q);
       end
   end
   C=sparse(ic, jc, vc, ra*rb, ca*cb);
end
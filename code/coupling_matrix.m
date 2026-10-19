function C = coupling_matrix(M, alpha, p, c0, beta, phi0, lambda)
%COUPLING_MATRIX Banded symmetric-Toeplitz mutual-coupling matrix, ICC Eq. (6).
%   c_m(alpha) = c0*exp(-beta*m*alpha*d0/lambda)*exp(j*psi_m),
%   psi_m = 2*pi*m*alpha*d0/lambda + phi0, for m = 1..p; c_0 = 1 on the
%   diagonal, and c_m = 0 for m > p (coupling range p, as Algorithm 1 assumes).
d0 = lambda/2;
c  = zeros(1, M);  c(1) = 1;
for m = 1:min(p, M-1)
    psi    = 2*pi*m*alpha*d0/lambda + phi0;
    c(m+1) = c0*exp(-beta*m*alpha*d0/lambda)*exp(1j*psi);
end
C = toeplitz(c);   % symmetric Toeplitz
end

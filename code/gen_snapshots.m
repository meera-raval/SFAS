function X = gen_snapshots(A, P, N, snr_db)
%GEN_SNAPSHOTS X = A*diag(sqrt(P))*S + noise, circular complex Gaussian
%   sources and AWGN. SNR is per source relative to unit source power.
K  = size(A, 2);
S  = (randn(K, N) + 1j*randn(K, N))/sqrt(2);
S  = diag(sqrt(P(:))) * S;
s2 = 10^(-snr_db/10);
Nn = sqrt(s2/2)*(randn(size(A,1), N) + 1j*randn(size(A,1), N));
X  = A*S + Nn;
end

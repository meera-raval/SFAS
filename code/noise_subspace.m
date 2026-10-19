function Un = noise_subspace(R, K)
%NOISE_SUBSPACE Eigenvectors of the M-K smallest eigenvalues of R.
R = (R + R')/2;
[U, L] = eig(R);
[~, idx] = sort(real(diag(L)), 'descend');
Un = U(:, idx(K+1:end));
end

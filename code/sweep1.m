% Stage-1 DOA accuracy vs. coupling strength c0 and coupling range p
rng(3); lambda=1; M=32; ac=0.1; N=500; SNR=10;
th=[-40 -20 10 30]; r=[30 300 1000 5000]; K=4; T=40;
fprintf('%5s %3s %6s %5s  %s\n','c0','p','|c1|','fail','Stage-1 RMSE per source (deg)  | max');
for p = [1 2 3]
 for c0 = [0 0.1 0.33]
  Cc = coupling_matrix(M, ac, p, c0, 2, pi/4, lambda);
  e = nan(T,K);
  for t=1:T
    X = gen_snapshots(Cc*esg_steering(th, r, ac, M, lambda), ones(1,K), N, SNR);
    t1 = sfas_stage1(X, K, p, ac, -90:0.05:90);
    if numel(t1)==K, e(t,:) = t1 - th; end
  end
  ok = ~isnan(e(:,1)); rm = sqrt(mean(e(ok,:).^2));
  fprintf('%5.2f %3d %6.3f %5d  %s | %5.2f\n', c0, p, abs(Cc(1,2)), sum(~ok), sprintf('%6.2f',rm), max(rm));
 end
end
% Noiseless far-field-approximation bias of the 30 lambda source vs range
fprintf('\nNoiseless Stage-1 bias, single source at -40 deg, no coupling:\n');
for rr = [15 30 60 100 300 1000]
  X = esg_steering(-40, rr, ac, M, lambda)*ones(1,500) + 1e-4*(randn(M,500)+1j*randn(M,500));
  t1 = sfas_stage1(X, 1, 3, ac, -90:0.01:90);
  fprintf('  r = %5d lambda: bias = %6.3f deg\n', rr, t1+40);
end

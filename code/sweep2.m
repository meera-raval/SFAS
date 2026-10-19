rng(4); lambda=1; M=32; N=500; SNR=10; th=[-40 -20 10 30]; r=[30 300 1000 5000]; K=4; T=40; p=1;
fprintf('%6s %6s %5s  %s\n','alphac','|c1|','fail','Stage-1 RMSE (deg)        | 30-lambda noiseless bias');
for ac = [0.1 0.2 0.3 0.45]
  Cc = coupling_matrix(M, ac, p, 0.33, 2, pi/4, lambda); e = nan(T,K);
  for t=1:T
    t1 = sfas_stage1(gen_snapshots(Cc*esg_steering(th,r,ac,M,lambda),ones(1,K),N,SNR), K, p, ac, -90:0.05:90);
    if numel(t1)==K, e(t,:)=t1-th; end
  end
  X = esg_steering(-40,30,ac,M,lambda)*ones(1,500)+1e-4*(randn(M,500)+1j*randn(M,500));
  b = sfas_stage1(X,1,p,ac,-90:0.01:90)+40; ok=~isnan(e(:,1));
  fprintf('%6.2f %6.3f %5d  %s | %6.2f\n', ac, abs(Cc(1,2)), sum(~ok), sprintf('%6.2f',sqrt(mean(e(ok,:).^2))), b);
end

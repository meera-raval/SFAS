% cfg: 'B' = ICC/journal stated spacing (ac=0.1, ae=1) + journal grids
%      'C' = spacing from journal figures (dc=0.1, de=1 -> ac=0.2, ae=2) + journal grids
lambda=1; M=32; N=500; SNR=10; th=[-40 -20 10 30]; r=[30 300 1000 5000]; K=4; p=1;
if cfg=='B', ac=0.1; ae=1; else, ac=0.2; ae=2; end
c0 = 0.3/exp(-2*ac*0.5);                     % |c1(ac)| = 0.30 (ICC), beta = 2
Cc = coupling_matrix(M, ac, p, c0, 2, pi/4, lambda);
opt = struct('r_grid', logspace(log10(5),4,216), 'dtheta_win',5, 'dr_win',60, ...
   'dtheta_step',0.015, 'dr_step',0.1, 'normalize',true, 'init','paper', 'n_pass',1);
fprintf('cfg %s: ac=%.1f ae=%.0f |c1(ac)|=%.3f |c1(ae)|=%.3f  Fraunhofer(ext)=%.0f lambda\n', cfg, ac, ae, ...
   abs(Cc(1,2)), c0*exp(-2*ae*0.5), 2*((M-1)*ae*0.5)^2);
rng(10);
for t = 1:T
  t1 = sfas_stage1(gen_snapshots(Cc*esg_steering(th,r,ac,M,lambda),ones(1,K),N,SNR), K, p, ac, -90:0.1:90);
  if numel(t1)<K, fprintf('trial %d: stage-1 unresolved\n',t); continue; end
  est = sfas_stage2(gen_snapshots(esg_steering(th,r,ae,M,lambda),ones(1,K),N,SNR), K, t1, ae, lambda, opt);
  fprintf('trial %2d | th1 err %s | th2 err %s | r err%% %s\n', t, sprintf('%6.2f',t1-th), ...
     sprintf('%7.3f',est.theta-th), sprintf('%7.1f',100*(est.r-r)./r));
end

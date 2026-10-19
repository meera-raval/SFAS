%SFAS_ICC_MC  Monte-Carlo check of the single-shot scenario (not in the ICC
%   paper, which shows only the single-shot example): RMSE over trials.
if ~exist('init_mode','var'), init_mode = 'paper'; end
rng(2);
lambda = 1; M = 32; alpha_c = 0.1; alpha_e = 1.0; N = 500; SNR = 10;
theta = [-40 -20 10 30]; r = [30 300 1000 5000]; K = 4; P = ones(1,K);
p = 3; Cc = coupling_matrix(M, alpha_c, p, 0.33, 2, pi/4, lambda);
opt.r_grid = logspace(log10(5), 4, 3000); opt.dtheta_win = 5;
opt.dr_win = 0.2*1e4; opt.dtheta_step = 0.01; opt.nr_loc = 800; opt.normalize = true; opt.init = init_mode; if ~exist('n_pass','var'), n_pass = 1; end; opt.n_pass = n_pass;
if ~exist('T','var'), T = 100; end; E1 = zeros(T,K); E2 = zeros(T,K); Er = zeros(T,K);
for t = 1:T
    Xc = gen_snapshots(Cc*esg_steering(theta, r, alpha_c, M, lambda), P, N, SNR);
    th1 = sfas_stage1(Xc, K, p, alpha_c, -90:0.1:90);
    if numel(th1) < K, E1(t,:)=NaN; E2(t,:)=NaN; Er(t,:)=NaN; continue; end
    Xe = gen_snapshots(esg_steering(theta, r, alpha_e, M, lambda), P, N, SNR);
    est = sfas_stage2(Xe, K, th1, alpha_e, lambda, opt);
    E1(t,:) = th1 - theta;  E2(t,:) = est.theta - theta;  Er(t,:) = (est.r - r)./r;
end
ok = ~isnan(E1(:,1)); fprintf('Mode: %s, passes %d | trials: %d | Stage-1 resolution failures: %d\n', init_mode, n_pass, T, sum(~ok));
E1=E1(ok,:); E2=E2(ok,:); Er=Er(ok,:);
fprintf('Stage-1 DOA RMSE (deg):  %s\n', sprintf('%8.3f', sqrt(mean(E1.^2))));
fprintf('Stage-2 DOA RMSE (deg):  %s\n', sprintf('%8.4f', sqrt(mean(E2.^2))));
fprintf('Range  median |err| (%%): %s\n', sprintf('%8.2f', 100*median(abs(Er))));
fprintf('Range  P(|err|<5%%):     %s\n', sprintf('%8.2f', mean(abs(Er)<0.05)));
fprintf('Coupling |c1| compressed = %.3f\n', abs(Cc(1,2)));

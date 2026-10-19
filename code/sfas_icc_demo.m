%SFAS_ICC_DEMO  Reproduces the single-shot mixed-field example (Sec. IV-A,
%   Figs. 2-3) of: T. Wu et al., "Scalable Fluid Antenna Systems for
%   Mixed-Field Source Localization," IEEE ICC 2026.
if ~exist('init_mode','var'), clear; init_mode = 'paper'; end
close all; rng(1);

%% Parameters (Sec. II-A, IV)
lambda  = 1;                 % work in wavelengths (paper: lambda = 1 mm)
M       = 32;
alpha_c = 0.1;  alpha_e = 1.0;
N       = 500;  SNR = 10;    % dB
theta   = [-40 -20 10 30];   % deg
r       = [30 300 1000 5000]*lambda;
K       = numel(theta);
P       = ones(1, K);

% Mutual coupling (Eq. (6)). Not fully specified in the paper: beta = 2 as
% stated, c0 chosen so |c1(alpha_c)| ~ 0.3 as stated, phi0 and p assumed.
p = 3;  beta = 2;  c0 = 0.33;  phi0 = pi/4;
Cc = coupling_matrix(M, alpha_c, p, c0, beta, phi0, lambda);

%% Stage 1: compressed configuration (true channel = coupling x exact ESG)
Xc = gen_snapshots(Cc*esg_steering(theta, r, alpha_c, M, lambda), P, N, SNR);
grid1 = -90:0.1:90;
[theta1, P1] = sfas_stage1(Xc, K, p, alpha_c, grid1);

%% Stage 2: extended configuration (coupling negligible, Eq. (8))
Xe = gen_snapshots(esg_steering(theta, r, alpha_e, M, lambda), P, N, SNR);
opt.r_grid      = logspace(log10(5), log10(1e4), 3000)*lambda;  % [rmin, rmax]
opt.dtheta_win  = 5;             % Delta_theta = 5 deg  (Sec. III-B)
opt.dr_win      = 0.2*max(opt.r_grid);   % Delta_r = 0.2 r_max  (Sec. III-B)
opt.dtheta_step = 0.01;
opt.nr_loc      = 800;
opt.normalize   = true;          % normalise ||a||; see notes
if ~exist('init_mode','var'), init_mode = 'paper'; end
opt.init        = init_mode;     % 'paper' or 'coarse2d'
if ~exist('n_pass','var'), n_pass = 1; end
opt.n_pass      = n_pass;        % 1 = paper; 2 = re-centre range window
est = sfas_stage2(Xe, K, theta1, alpha_e, lambda, opt);

%% Report
De = (M-1)*alpha_e*lambda/2;
fprintf('Extended aperture De = %.2f lambda, Fraunhofer 2De^2/lambda = %.0f lambda\n', De, 2*De^2/lambda);
fprintf('%4s %9s %9s %9s %10s %10s %10s %9s\n', 'k', 'theta', 'th1(S1)', 'th2(S2)', 'r_true', 'r_init', 'r_hat', 'r_err%');
for k = 1:numel(theta1)
    fprintf('%4d %9.2f %9.2f %9.3f %10.1f %10.1f %10.1f %9.2f\n', k, theta(k), theta1(k), ...
        est.theta(k), r(k), est.r_init(k), est.r(k), 100*abs(est.r(k)-r(k))/r(k));
end

%% Figures
figure; plot(grid1, 10*log10(P1/max(P1))); grid on; hold on;
xline_ = @(x) plot([x x], ylim, 'r--'); for k = 1:K, xline_(theta(k)); end
xlabel('\theta (deg)'); ylabel('Normalized spectrum (dB)'); title('Stage 1: compressed-array MUSIC');

figure;   % Fig. 2
for k = 1:K
    subplot(2, 2, k);
    semilogx(opt.r_grid, 10*log10(est.Pr{k}/max(est.Pr{k}))); hold on; grid on;
    yl = ylim; plot([r(k) r(k)], yl, 'r--');
    xlabel('Normalized distance (\lambda)'); ylabel('Spatial spectrum (dB)');
    title(sprintf('\\theta_%d^{(1)} = %.1f^o, r_%d = %g\\lambda', k, theta1(k), k, r(k)));
end

figure;   % Fig. 3
for k = 1:K
    subplot(2, 2, k);
    P2 = est.P2{k};
    surf(est.r_loc{k}, est.th_loc{k}, P2/max(P2(:)), 'EdgeColor', 'none');
    view(2); axis tight; colorbar;
    xlabel('Normalized distance (\lambda)'); ylabel('\theta (deg)');
    title(sprintf('Source %d: (%g^o, %g\\lambda)', k, theta(k), r(k)));
end

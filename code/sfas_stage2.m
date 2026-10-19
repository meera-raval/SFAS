function out = sfas_stage2(Xe, K, theta1, alpha_e, lambda, opt)
%SFAS_STAGE2 Joint DOA/range refinement in the extended configuration
%   (ICC paper Algorithm 2, Eqs. (12)-(15)).
%   opt fields:
%     r_grid      1-D range grid [rmin, rmax] for Eq. (14)
%     dtheta_win  Delta_theta (deg), dr_win Delta_r, dtheta_step (deg)
%     nr_loc      number of points of the local range grid (uniform in 1/r)
%     normalize   normalise ||a|| in the MUSIC denominator (true/false)
%     init        'paper'   : 1-D range search at fixed theta1 (Eq. (13)-(14))
%                 'coarse2d': coarse 2-D scan over |theta-theta1|<=Delta_theta
%                             and the full r_grid (robust initialiser, NOT in paper)
%     n_pass      number of (1-D range search + local 2-D) passes; 1 = paper.
%                 Pass >= 2 redoes Eq. (14) at the refined angle (NOT in paper)
[M, N] = size(Xe);
if ~isfield(opt, 'init'), opt.init = 'paper'; end
if ~isfield(opt, 'n_pass'), opt.n_pass = 1; end
Re = (Xe*Xe')/N;                                       % Eq. (12)
Un = noise_subspace(Re, K);
if opt.normalize
    spec = @(A) 1 ./ real(sum(abs(Un'*A).^2, 1) ./ sum(abs(A).^2, 1));
else
    spec = @(A) 1 ./ real(sum(abs(Un'*A).^2, 1));      % Eqs. (13),(15) as written
end
rg = opt.r_grid(:).';
Kh = numel(theta1);
out.theta = zeros(1, Kh); out.r = zeros(1, Kh); out.r_init = zeros(1, Kh); out.th_init = zeros(1, Kh);
for k = 1:Kh
    % Step 1: initial range, Eqs. (13)-(14)
    Pr = spec(esg_steering(theta1(k), rg, alpha_e, M, lambda));
    switch opt.init
        case 'paper'
            [~, i0] = max(Pr);  r0 = rg(i0);  t0 = theta1(k);
        case 'coarse2d'
            thc = theta1(k) + (-opt.dtheta_win:0.25:opt.dtheta_win);
            rc  = rg(1:5:end);
            Pc  = zeros(numel(thc), numel(rc));
            for it = 1:numel(thc), Pc(it,:) = spec(esg_steering(thc(it), rc, alpha_e, M, lambda)); end
            [~, idx] = max(Pc(:)); [it, ir] = ind2sub(size(Pc), idx);
            r0 = rc(ir);  t0 = thc(it);
    end
    for pass = 1:opt.n_pass
        if pass > 1   % re-run Eq. (14) at the refined angle
            Pr = spec(esg_steering(t0, rg, alpha_e, M, lambda));
            [~, i0] = max(Pr);  r0 = rg(i0);
        end
        % Step 2: localized 2-D refinement, Eq. (15)
        th_loc = t0 + (-opt.dtheta_win:opt.dtheta_step:opt.dtheta_win);
        rlo = max(rg(1), r0 - opt.dr_win);  rhi = min(rg(end), r0 + opt.dr_win);
        if isfield(opt, 'dr_step')        % linear grid, journal: delta_r = 0.1 lambda
            r_loc = rlo:opt.dr_step:rhi;
        else
            r_loc = 1 ./ linspace(1/rhi, 1/rlo, opt.nr_loc);   % uniform in 1/r
        end
        P2 = zeros(numel(th_loc), numel(r_loc));
        for it = 1:numel(th_loc)
            P2(it, :) = spec(esg_steering(th_loc(it), r_loc, alpha_e, M, lambda));
        end
        [~, idx] = max(P2(:));  [it, ir] = ind2sub(size(P2), idx);
        if pass < opt.n_pass, t0 = th_loc(it); end
    end
    out.theta(k) = th_loc(it);  out.r(k) = r_loc(ir);
    out.r_init(k) = r0;  out.th_init(k) = t0;
    out.Pr{k} = Pr;  out.P2{k} = P2;  out.th_loc{k} = th_loc;  out.r_loc{k} = r_loc;
end
end

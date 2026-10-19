function [theta1, P1, grid] = sfas_stage1(Xc, K, p, alpha_c, grid)
%SFAS_STAGE1 Coarse DOA estimation in the compressed configuration
%   (ICC paper Algorithm 1, Eqs. (9)-(11)).
%   Xc : M x N data from the compressed array (with mutual coupling)
%   p  : coupling range; the central M-2p elements are kept
[M, N] = size(Xc);
F   = [zeros(M-2*p, p), eye(M-2*p), zeros(M-2*p, p)];  % Eq. (9)
Xt  = F*Xc;                                            % step 2
Rt  = (Xt*Xt')/N;                                      % step 3
Un  = noise_subspace(Rt, K);                           % step 4
At  = ff_steering(grid, alpha_c, M);
At  = At(p+1:M-p, :);                                  % centre M-2p rows
P1  = 1 ./ real(sum(abs(Un'*At).^2, 1));               % Eq. (11)
loc = pick_peaks(P1, K);                               % step 8
theta1 = grid(loc);
end

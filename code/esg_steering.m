function A = esg_steering(theta_deg, r, alpha, M, lambda)
%ESG_STEERING Exact-spatial-geometry steering vectors, ICC paper Eqs. (2),(4),(5).
%   A = esg_steering(theta_deg, r, alpha, M, lambda)
%   theta_deg, r : scalars or equal-length vectors (one column per (theta,r) pair),
%                  or a scalar theta with a vector r (and vice versa).
%   Returns M x L matrix. Element 1 is the phase/range reference (p_1 = 0),
%   exactly as Eqs. (2) and (4) are written.
d0  = lambda/2;
p   = (0:M-1).' * alpha * d0;                 % Eq. (2)
th  = theta_deg(:).' * pi/180;
r   = r(:).';
if isscalar(th), th = th*ones(size(r)); end
if isscalar(r),  r  = r*ones(size(th)); end
rm  = sqrt(r.^2 + p.^2 - 2*p.*(r.*sin(th)));  % Eq. (4)  (M x L, implicit expansion)
phi = 2*pi/lambda * (rm - r);                 % excess path
A   = (r ./ rm) .* exp(1j*phi);               % Eq. (5)
end

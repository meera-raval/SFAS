function A = ff_steering(theta_deg, alpha, M)
%FF_STEERING Far-field steering vectors a_FF(theta, alpha) of Sec. II-C.
%   [1, e^{-j pi alpha sin(theta)}, ..., e^{-j (M-1) pi alpha sin(theta)}]^T
m = (0:M-1).';
A = exp(-1j*pi*alpha*m*sin(theta_deg(:).'*pi/180));
end

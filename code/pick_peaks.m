function [locs, vals] = pick_peaks(y, K)
%PICK_PEAKS Indices of the K largest local maxima of vector y (no toolboxes).
y  = y(:);
n  = numel(y);
ip = find([y(1) > y(2); y(2:n-1) >= y(1:n-2) & y(2:n-1) > y(3:n); y(n) > y(n-1)]);
[~, o] = sort(y(ip), 'descend');
ip = ip(o(1:min(K, numel(o))));
locs = sort(ip);
vals = y(locs);
end

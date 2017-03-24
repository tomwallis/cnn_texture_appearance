function C = crosscorr(a, b)

if nargin == 1
    a = (a - mean(a(:))) / std(a(:), 1);
    b = a;
else
    a = (a - mean(b(:))) / std(b(:), 1);
    b = (b - mean(b(:))) / std(b(:), 1);
end
[Ma, Na] = size(a);
[Mb, Nb] = size(b);
m = [1:Mb, repmat(Mb, 1, Ma-Mb-1), Mb:-1:1];
n = [1:Nb, repmat(Nb, 1, Na-Nb-1), Nb:-1:1];
N = m' * n;
C = xcorr2(a, b) ./ N;

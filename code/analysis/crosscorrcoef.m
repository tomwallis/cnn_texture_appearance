function R = crosscorrcoef(a, b)

[m, n] = size(b);
N = numel(b);
box = ones([m, n]) / N;
Ea = conv2(a, box, 'valid');
Eaa = conv2(a .* a, box, 'valid');
SDa = sqrt(Eaa - (Ea .* Ea));
Eb = mean(b(:));
SDb = std(b(:));
valid = @(x) x(m : end-m+1, n : end-n+1);
C = valid(xcorr2(a, b)) / N;
R = (C - Ea * Eb) ./ (SDa * SDb);

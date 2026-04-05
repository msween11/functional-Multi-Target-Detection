function [recspace, recfft] = kot_multipath(params, grids, psihat, psireg, ps, varargin)

use_simps = ismember('simps', varargin);
use_mixed = ismember('mixed', varargin);

mid = length(grids.Wext)/2;
Wpos = grids.Wext(mid+1:end);

empx = d14o(psihat, grids.dw, 1);

bispec_ratio = empx./psireg; %must use psihat here, psireg does not work

integrand_reg = diag(bispec_ratio)';

f_reg = integrand_reg(mid+1:end);

rec_reg = cumsimps(Wpos, f_reg);

g = zeros(1, length(Wpos));
g(1) = 0;
if use_mixed
    g(2) = imag(rec_reg(2));
else
    g(2) = rec_reg(2);
end

for k = 3:mid
    estimates = zeros(1, k);

    if use_mixed
        estimates(1) = imag(rec_reg(k));
    else
        estimates(1) = rec_reg(k);
    end

    for c = 1:k-1
        col  = mid+c;
        rows = mid+c:mid+k;
        if use_simps
            freq_slice = Wpos(c:k);
            integral = simps(freq_slice, bispec_ratio(rows, col));
        else
            integral = grids.dw * sum(bispec_ratio(rows, col));
        end
        if use_mixed
            estimates(c+1) = g(c) + imag(integral) + g(k-c);
        else
            estimates(c+1) = g(c) + integral - conj(g(k-c));
        end
    end

    g(k) = mean(estimates);
end

if use_mixed
    recmp = sqrt(ps(mid+1:end)) .* exp(1i * g);
else
    recmp = exp(g);
end


recfft = [conj(flip(recmp)) (recmp)];


%hfix recovery in space now
phiK = @(w) 1 .* (-1 <= w & w <= 1);


rec = zeros(1, length(grids.Dext));
for j = 1:length(grids.Dext) 
    rec(j) = (1/(2*pi)) * simps( ...
    grids.Wext, ...
    exp(1i*grids.Wext*grids.Dext(j)) .* ...
    recfft .* phiK(params.h*grids.Wext));
end
t = length(grids.D)/2;  
tt = length(grids.Dext)/2;
rec_trunc = real(rec(tt-t+1:tt+t)); 
recspace = rec_trunc/max(abs(rec_trunc));


end

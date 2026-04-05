function [recspace, recfft] = kotlarski_auto(params, grids, psihat, psireg, varargin)
%KOTLARSKI_AUTO  Unified Kotlarski characteristic-function recovery.
%
%   [recspace, recfft] = kotlarski_auto(params, grids, psihat, psireg)
%   [recspace, recfft] = kotlarski_auto(..., 'method',    METHOD)
%   [recspace, recfft] = kotlarski_auto(..., 'ps',        ps)
%   [recspace, recfft] = kotlarski_auto(..., 'multipath', true)
%   [recspace, recfft] = kotlarski_auto(..., 'use_reg',   true)
%   [recspace, recfft] = kotlarski_auto(..., 'gradx',     psihat_gradx)
%
%   Required inputs
%     params   - struct with fields: zero_cond, thresh, h
%     grids    - struct with fields: Wext, dw, Dext, D
%     psihat   - empirical joint characteristic function matrix
%     psireg   - regularised joint characteristic function matrix
%
%   Optional name-value pairs
%     'method'    - quadrature rule for all integration steps:
%                     'trapz'    (default) cumtrapz / trapz
%                     'simpsons' cumsimps  / simps
%                     'riemann'  left Riemann sum (matches original kotlarski.m)
%     'ps'        - spectral density vector; enables mixed mode: amplitude
%                   taken from sqrt(ps), only phase is integrated
%     'multipath' - logical (default false); averages estimates across all
%                   bispectrum columns at each frequency step
%     'use_reg'   - logical (default false); use psireg as the integrand
%                   denominator instead of psihat
%     'gradx'     - precomputed horizontal gradient of psihat (4th output of
%                   get_psihat_auto); replaces d14o(psihat) with the exact
%                   FT-based gradient FT[+i*z1*A3], which is more accurate
%                   at high noise

p = inputParser;
addParameter(p, 'ps',        []);
addParameter(p, 'method',    'trapz', @(x) ismember(x, {'trapz','simpsons','riemann'}));
addParameter(p, 'multipath', false);
addParameter(p, 'use_reg',   true);
addParameter(p, 'gradx',     []);
parse(p, varargin{:});

ps        = p.Results.ps;
method    = p.Results.method;
multipath = logical(p.Results.multipath);
use_mixed = ~isempty(ps);
use_reg   = logical(p.Results.use_reg);

mid  = length(grids.Wext) / 2;
Wpos = grids.Wext(mid+1:end);

if ~isempty(p.Results.gradx)
    empx = p.Results.gradx;
else
    empx = d14o(psihat, grids.dw, 1);
end
adx  = diag(empx);

% Denominator for single-path integrand: psihat diagonal or psireg diagonal
if use_reg
    ad_denom = diag(psireg);
else
    ad_denom = diag(psihat);
end

%% ----------------------------------------------------------------
%  PHASE / AMPLITUDE RECOVERY
%% ----------------------------------------------------------------

if multipath
    %------------------------------------------------------------------
    % Multipath: at each frequency k, average estimates from all
    % bispectrum columns plus the diagonal path.
    %------------------------------------------------------------------
    if use_reg
        bispec_ratio = empx ./ psireg;
    else
        bispec_ratio = empx ./ psihat;
    end

    integrand_diag = diag(bispec_ratio)';
    f_diag         = integrand_diag(mid+1:end);
    rec_diag       = cumint(Wpos, f_diag, method, grids.dw);

    g    = zeros(1, mid);
    g(1) = 0;
    if use_mixed
        g(2) = imag(rec_diag(2));
    else
        g(2) = rec_diag(2);
    end

    for k = 3:mid
        estimates    = zeros(1, k);
        if use_mixed
            estimates(1) = imag(rec_diag(k));
        else
            estimates(1) = rec_diag(k);
        end

        for c = 1:k-1
            col      = mid + c;
            rows     = (mid+c):(mid+k);
            integral = int1d(Wpos(c:k), bispec_ratio(rows, col), method, grids.dw);
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
    recfft = [conj(flip(recmp)) recmp];

    if params.zero_cond == 1
        if use_mixed
            ps_zf = ps;
        else
            ps_zf = real(diag(psihat))';
        end
        rec_pos = zero_flip(recfft(mid+1:end), ps_zf, params, grids, mid);
        recfft  = [conj(flip(rec_pos)) rec_pos];
    end

elseif use_mixed
    %------------------------------------------------------------------
    % Mixed single-path: integrate only the phase, amplitude from sqrt(ps).
    %------------------------------------------------------------------
    integrand     = (adx ./ ad_denom)';
    phase_deriv   = imag(integrand(mid+1:end));
    phase_pos     = cumint(Wpos, phase_deriv, method, grids.dw);
    rec_pos       = ps(mid+1:end).^(0.5) .* exp(1i * phase_pos);

    if params.zero_cond == 1
        rec_pos = zero_flip(rec_pos, ps, params, grids, mid);
    end
    recfft = [conj(flip(rec_pos)) rec_pos];

else
    %------------------------------------------------------------------
    % Standard single-path: integrate full complex log-derivative.
    %------------------------------------------------------------------
    integrand = (adx ./ ad_denom)';
    f         = integrand(mid+1:end);
    rec_pos   = exp(cumint(Wpos, f, method, grids.dw));

    if params.zero_cond == 1
        ps_zf   = real(diag(psihat))';
        rec_pos = zero_flip(rec_pos, ps_zf, params, grids, mid);
    end
    recfft = [conj(flip(rec_pos)) rec_pos];
end

%% ----------------------------------------------------------------
%  RECOVERY IN SPACE  (inverse Fourier transform)
%% ----------------------------------------------------------------
phiK = @(w) 1 .* (-1 <= w & w <= 1);
rec  = zeros(1, length(grids.Dext));
for j = 1:length(grids.Dext)
    integrand_fft = exp(1i * grids.Wext * grids.Dext(j)) ...
                    .* recfft .* phiK(params.h * grids.Wext);
    rec(j) = (1/(2*pi)) * int1d(grids.Wext, integrand_fft, method, grids.dw);
end

t  = length(grids.D)    / 2;
tt = length(grids.Dext) / 2;
rec_trunc = real(rec(tt-t+1 : tt+t));
recspace  = rec_trunc / max(abs(rec_trunc));

end

%% ================================================================
%  LOCAL HELPERS
%% ================================================================

function rec_pos = zero_flip(rec_pos, ps_full, params, grids, mid)
%ZERO_FLIP  Detect spectral zeros in ps_full and flip sign of rec_pos
%           (positive-frequency half, length mid) across each zero
%           Also zeros out a small window around each detected zero
growth_filt = .001;
g1          = d14o(ps_full, grids.dw, 2);
g2          = d24o(ps_full, grids.dw, 2);
pidz        = findz(g1);
pidz_filt   = pidz(abs(g2(pidz)) > growth_filt);
pidz_filt   = pidz_filt(grids.Wext(pidz_filt) > 0);
pz_filt     = grids.Wext(pidz_filt);
fliploc     = pz_filt(abs(ps_full(pidz_filt)) < params.thresh);
[~, loc]    = ismember(fliploc, grids.Wext);
loc         = loc - mid;
loc         = [loc mid];

for j = 1:length(loc)-1
    rec_pos(loc(j)+1:loc(j+1)) = (-1)^j * rec_pos(loc(j)+1:loc(j+1));
end

ep = 2; zloc = [];
for w = 1:ep
    zloc = [zloc loc-w loc loc+w];
end
zloc = unique(zloc);
zloc(end-ep+1:end) = [];
zloc = zloc(zloc >= 1 & zloc <= mid);
rec_pos(zloc) = 0;
end

function z = cumint(x, y, method, dw)
%CUMINT  Cumulative integral of y with respect to x
switch method
    case 'trapz'
        z = cumtrapz(x, y);
    case 'simpsons'
        z = cumsimps(x, y);
    case 'riemann'
        % Matches original kotlarski.m: rec(1)=0, rec(k)=dw*sum(f(2:k))
        z = dw * [0, cumsum(y(2:end))];
end
end

function val = int1d(x, y, method, dw)
%INT1D  Definite integral of y with respect to x
switch method
    case 'trapz'
        val = trapz(x, y);
    case 'simpsons'
        val = simps(x, y);
    case 'riemann'
        val = dw * sum(y);
end
end

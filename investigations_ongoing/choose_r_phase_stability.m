function [r_opt, diag] = choose_r_phase_stability(psihat, grids, params, rvals)
%CHOOSE_R_PHASE_STABILITY  Choose regularization r via phase stability
%
%   r_opt = choose_r_phase_stability(psihat, grids, params)
%   [r_opt, diag] = choose_r_phase_stability(...)
%
% Inputs:
%   psihat  : empirical characteristic function (vector over grids.Wext)
%   grids   : struct with fields Wext, dw
%   params  : struct with fields n, h, L
%   rvals   : (optional) vector of candidate r values
%
% Output:
%   r_opt   : chosen regularization constant
%   diag    : diagnostic struct (phase variation curve, etc.)

    if nargin < 4 || isempty(rvals)
        rvals = logspace(-3, 0, 30);   % default scan
    end

    W = grids.Wext(:);
    dw = grids.dw;
    n  = params.n;
    L  = grids.L;
    h  = params.h;

    % frequency mask: ONLY where psi_h will be used
    mask = abs(W) <= 1/h;

    V = zeros(size(rvals));

    for k = 1:length(rvals)
        r = rvals(k);

        % regularize (same rule as your code)
        psireg = psihat ./ ...
            (L * min(ones(size(psihat))/L, r*sqrt(n)*abs(psihat)));

        % phase on kept frequencies only
        phi = unwrap(angle(psireg(mask)));

        % phase-variation diagnostic
        V(k) = sum(abs(diff(phi))) * dw;
    end

    % detect "elbow" via max curvature in log-space
    logV = log(V);
    d2 = diff(logV,2);
    [~, idx] = max(d2);

    % +1 because diff(logV,2) shifts indices
    r_opt = rvals(idx+1);

    % diagnostics
    diag.rvals = rvals;
    diag.phase_variation = V;
    diag.log_phase_variation = logV;
    diag.curvature = d2;
    diag.idx = idx+1;
end

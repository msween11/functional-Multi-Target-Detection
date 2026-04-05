function [recspace, recfft] = kotlarski_simps_mixed(params, grids, psihat, psireg, ps)

%% DOING KOTLARSKI NOW

mid = length(grids.Wext)/2;
Wpos = grids.Wext(mid+1:end);

empx = d14o(psihat, grids.dw, 1); 
adx = diag(empx); 
ad_reg = diag(psireg);
integrand_reg = adx./ad_reg; integrand_reg = integrand_reg';

%new additions for phase construction
phase_deriv_pos = imag(integrand_reg(mid+1:end));
phase_pos = cumsimps(Wpos, phase_deriv_pos);
phase = [-fliplr(phase_pos) phase_pos];

rec1_reg = ps.^(.5).* exp(1i * phase);

if(params.zero_cond == 1)
    growth_filt = .001;
    g  = d14o(ps, grids.dw, 2);
    gg = d24o(ps, grids.dw, 2);

    pidz = findz(g);
    pidz_filt = pidz(abs(gg(pidz)) > growth_filt);
    pidz_filt = pidz_filt(grids.Wext(pidz_filt) > 0);
    pz_filt   = grids.Wext(pidz_filt);
    fliploc   = pz_filt(abs(ps(pidz_filt)) < params.thresh);

    [~, loc] = ismember(fliploc, grids.Wext);
    loc = loc - mid;
    loc = [loc mid];

    % flipping positive and negative frequencies symmetrically
    for j = 1:length(loc)-1
        pos_idx = mid + loc(j)+1 : mid + loc(j+1);
        neg_idx = mid + 1 - loc(j+1) : mid - loc(j);
        rec1_reg(pos_idx) = (-1)^j * rec1_reg(pos_idx);
        rec1_reg(neg_idx) = (-1)^j * rec1_reg(neg_idx);
    end

    % zeroing windows around zeros, symmetrically
    ep = 2; zloc = [];
    for w = 1:ep
        temp = [loc-w loc loc+w];
        zloc = [zloc temp];
    end
    clear('temp');
    zloc = unique(zloc);
    zloc(end-ep+1:end) = [];
    zloc = zloc(zloc >= 1 & zloc <= mid);

    rec1_reg(mid + zloc) = 0;
    rec1_reg(mid + 1 - zloc) = 0;
end

recfft = rec1_reg;

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
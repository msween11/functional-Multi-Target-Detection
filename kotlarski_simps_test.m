function [recspace, recfft] = kotlarski_simps_test(params, grids, psihat)

%% DOING KOTLARSKI NOW

mid = length(grids.Wext)/2;
Wpos = grids.Wext(mid+1:end);

empx = d14o(psihat, grids.dw, 1); %this can be precomputed before each path

adx = diag(empx); %start for each path
ad_emp = diag(psihat); %indexing points from psihat and gradient psihat

integrand_emp = adx./ad_emp; integrand_emp = integrand_emp'; %basic computation

f_emp = integrand_emp(mid+1:end); %indexing
rec_emp = cumsimps(Wpos, f_emp); %only real computational line

psireg = psihat./(grids.L*min(ones(length(grids.Wext))/grids.L,...
    params.r*sqrt(params.n)*abs(psihat)));

ad_reg = diag(psireg);
integrand_reg = adx./ad_reg; integrand_reg = integrand_reg';

f_reg = integrand_reg(mid+1:end);
rec_reg = cumsimps(Wpos, f_reg);
rec1_reg = (exp(rec_reg));  %end for each path

if(params.zero_cond == 1)
    growth_filt = .001;
    ps = real(ad_emp)';
    g = d14o(ps, grids.dw,2);
    gg = d24o(ps, grids.dw,2);

    pidz = findz(g);
    pidz_filt = pidz(abs(gg(pidz))> growth_filt );
    pidz_filt = pidz_filt(grids.Wext(pidz_filt)>0);
    pz_filt = grids.Wext(pidz_filt);
    fliploc = pz_filt(abs(ps(pidz_filt))< params.thresh); %this is by thresh

    %getting the index
    [~,loc] = ismember(fliploc, grids.Wext);
    loc = loc-length(grids.Wext)/2;
    loc = [loc length(grids.Wext)/2];

    %flipping
    for j = 1:length(loc)-1
        rec1_reg(loc(j)+1:loc(j+1))=(-1)^j*rec1_reg(loc(j)+1:loc(j+1));
    end

    %zeroing the little windows around the zeros
    ep = 2; zloc = [];
    for w=1:ep
        temp = [loc-w loc loc+w];
        zloc = [zloc temp];
    end
    clear('temp');
    zloc = unique(zloc);
    zloc(end-ep+1:end) = [];     
    rec1_reg(zloc) = 0;
end


recfft = [conj(flip(rec1_reg)) (rec1_reg)];


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
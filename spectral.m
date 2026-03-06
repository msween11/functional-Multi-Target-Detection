function [spectral_space, spectral_freq] = spectral(params, grids, psireg,f)

dd = grids.lenD/2;
Dd = length(grids.Dext)/2;
p = circshift(real(diag(psireg)), Dd); 

phases = get_B_phase((psireg));
phases2 = get_phase_from_bispectrum_gap(phases, size(psireg,1));
spectral_freq = sqrt(p).*phases2;
spectral_freq = circshift(spectral_freq', Dd);

spectral_space = zeros(1, length(grids.Dext));

phiK = @(w) 1 .* (-1 <= w & w <= 1);

for j = 1:length(grids.Dext)
    spectral_space(j) = (grids.dw/(2*pi)).* sum(exp(1i.*grids.Wext*grids.Dext(j)).*...
        spectral_freq .* phiK(params.h*grids.Wext)); %simpsons rule here
end

spectral_space = real(spectral_space)./max(abs(real(spectral_space)));

spectral_space = align_to_reference(spectral_space, f(grids.Dext));
spectral_space = spectral_space(Dd-dd+1:Dd+dd)'; 

end
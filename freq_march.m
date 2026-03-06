function [FM_space, FM_freq] = freq_march(params, grids, psireg,f)

dd = grids.lenD/2;
Dd = length(grids.Dext)/2;
p = circshift(real(diag(psireg)), Dd); %choose psihat here maybe?

phases = phases_from_bispectrum_FM_real(center_BS(psireg));
FM_freq = sqrt(p).*phases;
FM_freq = circshift(FM_freq', Dd);

FM_space = zeros(1, length(grids.Dext));

phiK = @(w) 1 .* (-1 <= w & w <= 1);

for j = 1:length(grids.Dext)
    FM_space(j) = (grids.dw/(2*pi)).* sum(exp(1i.*grids.Wext*grids.Dext(j)).*...
        FM_freq .* phiK(params.h*grids.Wext)); %simpsons rule here
end
    
FM_space = real(FM_space)./max(abs(real(FM_space)));

FM_space = align_to_reference(FM_space, f(grids.Dext));
FM_space = FM_space(Dd-dd+1:Dd+dd)'; 

end
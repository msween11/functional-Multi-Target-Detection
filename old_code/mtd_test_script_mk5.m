rng(1);

[params, grids] = setparams(n= 2^10, l = 5, B = 10, z = 1, ...
    t = .01, sigma = .1, lambda = .1, r = .1, h = .01);

func_choice = 'f1';
f = choose_function(func_choice); 

params.zero_cond = strcmp(func_choice, 'f4');

datastrip = get_obs(params, grids, f);

  %avoid streamed for now!
[ps, psihat,psireg] = get_psihat_auto(params, grids, f, datastrip);
[rec_kot, recfft_kot] = kotlarski_simps_test(params, grids, psihat, psireg);
[rec_mixed, recfft_mixed] = kotlarski_simps_mixed(params,grids, psihat, psireg, ps);
[FM_space, FM_freq] = freq_march(params, grids, psireg,f);
[spectral_space, spectral_freq] = spectral(params, grids, psireg,f);

%% FINAL PLOTTING

recspace = align_to_reference(recspace, f(grids.D))';
plotting = symfft(f(grids.Dext));
plotting = plotting ./plotting(length(grids.Dext)/2+1);

% l8_freq_rel_err = max(abs(recfft - plotting))/max(abs(plotting));
l8_space_err = max(abs(recspace-f(grids.D)))/max(abs(f(grids.D)));
l8_fm_err = max(abs(FM_space-f(grids.D)))/max(abs(f(grids.D)));
l8_spec_err = max(abs(spectral_space-f(grids.D)))/max(abs(f(grids.D)));

% figure
% plot(grids.Wext,real(plotting))
% hold on
% plot(grids.Wext, real(recfft))
% 
% figure
% plot(grids.D, recspace)
% hold on
% plot(grids.D,FM_space)
% plot(grids.D,spectral_space) 
% plot(grids.D,(f(grids.D)))
% legend(sprintf('Kot. = %.2f', l8_space_err),...
%     sprintf('FM = %.2f', l8_fm_err),...
%     sprintf('Spec. = %.2f', l8_spec_err),...
%     'True')

mid = length(grids.Wext)/2;
Wpos = grids.Wext(mid+1:end);

empx = d14o(psihat, grids.dw, 1);

bispec_ratio = empx./psireg; %must use psihat here, psireg does not work

integrand_reg = diag(bispec_ratio)';

f_reg = integrand_reg(mid+1:end);

rec_reg = cumsimps(Wpos, f_reg);
rec1_reg = exp(rec_reg);

recfft = [conj(flip(rec1_reg)) (rec1_reg)];

g = zeros(1, length(Wpos));
g(1) = 0;
g(2) = rec_reg(2);   
%%
for k = 3:mid
    estimates = zeros(1, k);

    % Estimate 1: diagonal Kotlarski — k-th positive frequency
    estimates(1) = rec_reg(k);

    % Estimates 2 to k: with c = 1, ..., k-1
    for c = 1:k-1
        col  = mid+c;
        rows = mid+c:mid+k;
        freq_slice = Wpos(c:k);
        % integral = simps(freq_slice, bispec_ratio(rows, col));
        integral = grids.dw*sum(bispec_ratio(rows, col));
        estimates(c+1) = (g(c)) + conj(integral) - conj(g(k-c));
    end

    g(k) = mean(estimates);
end

recmp = exp(g);
%%
plotting = imag(symfft(f(grids.Dext)));
plotting = plotting./max(plotting);
plotting = plotting(mid+1:end);

figure
plot(imag(recmp(1:200)))
hold on
plot(imag(rec1_reg(1:200)))
plot(plotting(1:200))
legend('mp', 'og', 'true')
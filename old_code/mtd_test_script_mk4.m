%mk4 updates: split code into smaller functions, 
%...implemented params and grids to simplify function inputs
%...implemented 4th order derivs in kotlarski
%... changed shifts to well_sep

%try higher freq signal
%compare freq recovery for algorithms instead, FM worse
%also compare l2 recovery in space, freq
%l2, l8 rec for space and freq for all algs PLUS
%also fix n, sigma, 1plot vary l, other B, all algs err

%ANNA meeting notes
%try increasing lambda to increase corr to mess with FM
%use A2M to check power spectrum error vs diag(psihat)!
%^^^for all algs
%l = 1 to 7
%B = 1 to 20
% rng(1);

[params, grids] = setparams(n= 2^10, l = 5, B = 10, z = 1, ...
    t = .01, sigma = .5, lambda = .1, r = .0036, h = .05);

func_choice = 'f1';
f = choose_function(func_choice); 

if strcmp(func_choice, 'f4')
   params.zero_cond = 1;
else 
   params.zero_cond = 0;
end

datastrip = get_obs(params, grids, f);

[ A3M, psihat, psireg] = get_psihat(params, grids, f, datastrip);

[r_opt, diag] = choose_r_phase_stability(psihat, grids, params);

params.r = r_opt;

[recspace, recfft] = kotlarski_simps_test(params, grids, psihat);

[FM_space, FM_freq] = freq_march(params, grids, psireg,f);

%% FINAL PLOTTING


recspace = align_to_reference(recspace, f(grids.D))';
plotting = symfft(f(grids.Dext));
plotting = plotting ./plotting(length(grids.Dext)/2+1);

l8_freq_rel_err = max(abs(recfft - plotting))/max(abs(plotting));
l8_space_err = max(abs(recspace-f(grids.D)))/max(abs(f(grids.D)));
l8_alg_err = max(abs(FM_space-f(grids.D)))/max(abs(f(grids.D)));


figure
plot(grids.Wext(length(grids.Wext)/2:end), ...
    real(plotting(length(plotting)/2:end)))
hold on
plot(grids.Wext(length(grids.Wext)/2:end), ...
    real(recfft(length(recfft)/2:end)))
% plot(Wext(length(Wext)/2:end), imag(recfft(length(recfft)/2:end)))

figure
plot(grids.D, recspace)
hold on
plot(grids.D,FM_space)
hold on
plot(grids.D,(f(grids.D)))
hold on
% plot(grids.D,spectral(params,grids, psireg, f))
legend('ours', 'FM', 'true', 'spectral')

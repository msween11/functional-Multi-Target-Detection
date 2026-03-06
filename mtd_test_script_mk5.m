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
rng(1);

[params, grids] = setparams(n= 2^10, l = 5, B = 10, z = 1, ...
    t = .01, sigma = .5, lambda = .1, r = .1, h = .05);

func_choice = 'f1';
f = choose_function(func_choice); 

if strcmp(func_choice, 'f4')
   params.zero_cond = 1;
else 
   params.zero_cond = 0;
end

datastrip = get_obs(params, grids, f);

  %avoid streamed for now!
[A1M, A2M,A3M, psihat,psireg] = get_psihat_auto(params, grids, f, datastrip);
n1 = grids.lenD;
n2 = length(grids.Dext);
nside = (n2-n1)/2;
A2M_padded = [zeros(1,nside) A2M zeros(1,nside)];
ps = real(symfft(A2M_padded));
ps = max(ps,0); ps = ps./max(ps); ps = sqrt(ps);

[recspace, recfft] = kotlarski_simps_test(params, grids, psihat);

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

figure
plot(grids.Wext,real(plotting))
hold on
plot(grids.Wext, real(recfft))

figure
plot(grids.D, recspace)
hold on
plot(grids.D,FM_space)
plot(grids.D,spectral_space) 
plot(grids.D,(f(grids.D)))
legend(sprintf('Kot. = %.2f', l8_space_err),...
    sprintf('FM = %.2f', l8_fm_err),...
    sprintf('Spec. = %.2f', l8_spec_err),...
    'True')

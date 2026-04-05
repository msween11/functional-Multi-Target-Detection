%mk6 updates: 
% main goal of script: 1) get est ps from A2M working AND
rng(1);

[params, grids] = setparams(n= 2^10, l = 5, B = 10, z = 1, ...
    t = .01, sigma = .0000001, lambda = .1, r = .1, h = .05);

covr = @(x) (params.sigma^2)*exp(-abs(x).^2./(2*params.lambda^2)); 
covr_freq = real(grids.dx*symfft(covr(grids.Dext)));
%both covr and covr_freq are scaled correctly according to cont theory

func_choice = 'f2';
f = choose_function(func_choice); 

if strcmp(func_choice, 'f4')
   params.zero_cond = 1;
else 
   params.zero_cond = 0;
end

% f = @(x) zeros(1, length(x));

datastrip = get_obs(params, grids, f, type='unif');

%avoid streamed for now!
[ps, psihat,psireg] = get_psihat_auto(params, grids, f, datastrip);
n1 = grids.lenD;
n2 = length(grids.Dext);
nside = (n2-n1)/2;

[recspace, recfft] = kotlarski_simps_test(params, grids, psihat, psireg);

[FM_space, FM_freq] = freq_march(params, grids, psireg,f);

[spectral_space, spectral_freq] = spectral(params, grids, psireg,f);

%% Testing mag / phase separated recovery

%true power spectrum: |\hatf(w)|^2
func = f(grids.Dext);
true_ps = abs(grids.dx*symfft(func)).^2;

%kot_ps is normalized estimated from abs of diag of psihat
kot_ps = abs(diag(psihat))';

%computing kot_ps and ps errors
ps_err = max(abs(ps - true_ps))/max(abs(true_ps));
kot_err = max(abs(kot_ps - true_ps))/max(abs(true_ps));

%plotting power spectrums
figure; hold on;
plot(true_ps)
plot(ps)
plot(kot_ps)
legend('true',sprintf('test, err = %f', ps_err), sprintf('kot, err = %f', kot_err))
%% reconstructing recfft from ps

%true fft CAUTION REAL() HERE !!!!
true_fft = real(grids.dx*symfft(f(grids.Dext)));

%phase information from third order information
mid = length(grids.Wext)/2;
Wpos = grids.Wext(mid+1:end);
empx = d14o(psihat, grids.dw, 1); %angle(psihat) here?
adx = diag(empx); 
ad_emp = diag(psihat);
integrand_emp = adx./ad_emp; integrand_emp = integrand_emp';
f_emp = integrand_emp(mid+1:end);
rec_emp = cumsimps(Wpos, f_emp);
ad_reg = diag(psireg);
integrand_reg = adx./ad_reg; integrand_reg = integrand_reg';

phase_deriv_pos = imag(integrand_reg(mid+1:end));
phase_pos = cumsimps(Wpos, phase_deriv_pos);
phase = [-fliplr(phase_pos) phase_pos];

phase2 = angle(recfft); 


%reconstructions from 3rd and 2nd order information
%kot_recon and recfft should match VERY closely
kot_recon = kot_ps.^(.5).* exp(1i * phase);
ps_recon = ps.^(.5).* exp(1i * phase);

%reconstruction errors
third_err = max(abs(recfft./max(abs(recfft)) - true_fft))/max(abs(true_fft));
ps_err = max(abs(ps_recon - true_fft))/max(abs(true_fft));
kot_err = max(abs(kot_recon - true_fft))/max(abs(true_fft));

figure; hold on;
plot(true_fft)
plot(real(recfft)./max(abs(recfft)))
plot(real(kot_recon));
plot(real(ps_recon));

legend('true ft',...
    sprintf('all 3rd order, kot. err = %f', third_err),...
    sprintf('3rd order recon. err = %f',kot_err),...
    sprintf('mixed 2nd 3rd order recon. err = %f',ps_err))

%%




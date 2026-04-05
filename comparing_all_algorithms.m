rng(1);

[params, grids] = setparams(n= 2^10, l = 5, B = 10, z = 1, ...
    t = .01, sigma = .01, lambda = .1, r = .01, h = .05);

func_choice = 'f5';
f = choose_function(func_choice);

params.zero_cond = strcmp(func_choice, 'f4');

datastrip = get_obs(params, grids, f);

%avoid streamed for now!
[ps, psihat, psireg, psihat_gradx] = get_psihat_auto(params, grids, f, datastrip);
[rec,       recfft]       = kotlarski_auto(params, grids, psihat, psireg, 'method', 'riemann', 'gradx', psihat_gradx);
[rec_simps, recfft_simps] = kotlarski_auto(params, grids, psihat, psireg, 'method', 'trapz',   'gradx', psihat_gradx);
[rec_mixed, recfft_mixed] = kotlarski_auto(params, grids, psihat, psireg, 'method', 'trapz',   'ps', ps, 'gradx', psihat_gradx);
[mp,        mp_fft]       = kotlarski_auto(params, grids, psihat, psireg, 'method', 'riemann',  'multipath', true, 'gradx', psihat_gradx);
[mp_sm,     mp_sm_fft]    = kotlarski_auto(params, grids, psihat, psireg, 'method', 'trapz',    'multipath', true, 'ps', ps, 'gradx', psihat_gradx);
[FM_space,     FM_freq]         = freq_march(params, grids, psireg, f);
[FM_mix_space, FM_mix_freq]    = freq_march(params, grids, psireg, f, ps);
[spec_space,   spec_freq]       = spectral(params, grids, psireg, f);

%% Shared reference and normalisation
mid   = length(grids.Wext)/2;
dc    = mid + 1;

true_fft_raw  = symfft(f(grids.Dext));
true_fft_norm = abs(true_fft_raw) / abs(true_fft_raw(dc));

norm_fft = @(x) abs(x) / abs(x(dc));
recfft_n    = norm_fft(recfft);
recfft_sn   = norm_fft(recfft_simps);
recfft_mn   = norm_fft(recfft_mixed);
mp_fft_n    = norm_fft(mp_fft);
mp_sm_fft_n = norm_fft(mp_sm_fft);
FM_fft_n     = norm_fft(FM_freq);
FM_mix_fft_n = norm_fft(FM_mix_freq);
spec_fft_n   = norm_fft(spec_freq);

rel_linf = @(x, ref) max(abs(x - ref)) / max(abs(ref));

%% Plot 1: Full-spectrum fft overlay
e1 = rel_linf(recfft_n,    true_fft_norm);
e2 = rel_linf(recfft_sn,   true_fft_norm);
e3 = rel_linf(recfft_mn,   true_fft_norm);
e4 = rel_linf(mp_fft_n,    true_fft_norm);
e5 = rel_linf(mp_sm_fft_n, true_fft_norm);
e6 = rel_linf(FM_fft_n,      true_fft_norm);
e7 = rel_linf(FM_mix_fft_n,  true_fft_norm);
e8 = rel_linf(spec_fft_n,    true_fft_norm);

figure;
plot(grids.Wext, true_fft_norm,  'k',  'LineWidth', 1.5); hold on;
plot(grids.Wext, recfft_n);
plot(grids.Wext, recfft_sn);
plot(grids.Wext, recfft_mn);
plot(grids.Wext, mp_fft_n);
plot(grids.Wext, mp_sm_fft_n);
plot(grids.Wext, FM_fft_n);
plot(grids.Wext, FM_mix_fft_n);
plot(grids.Wext, spec_fft_n);
legend('True', ...
    sprintf('Kot. (%.2f)',          e1), ...
    sprintf('Kot. simps (%.2f)',    e2), ...
    sprintf('Mixed (%.2f)',         e3), ...
    sprintf('Multipath (%.2f)',     e4), ...
    sprintf('MP simps+mix (%.2f)', e5), ...
    sprintf('FM (%.2f)',            e6), ...
    sprintf('FM mixed (%.2f)',      e7), ...
    sprintf('Spectral (%.2f)',      e8));
title('FFT Recovery — Full Spectrum');
xlabel('\omega'); ylabel('|F(w)| / |F(0)|');

%% Plot 2: First 200 nonneg frequencies
idx   = dc : dc + 199;
W200  = grids.Wext(idx);
true200 = true_fft_norm(idx);

e1_200 = rel_linf(recfft_n(idx),    true200);
e2_200 = rel_linf(recfft_sn(idx),   true200);
e3_200 = rel_linf(recfft_mn(idx),   true200);
e4_200 = rel_linf(mp_fft_n(idx),    true200);
e5_200 = rel_linf(mp_sm_fft_n(idx), true200);
e6_200 = rel_linf(FM_fft_n(idx),      true200);
e7_200 = rel_linf(FM_mix_fft_n(idx), true200);
e8_200 = rel_linf(spec_fft_n(idx),   true200);


figure;
plot(W200, true200,            'k',  'LineWidth', 1.5); hold on;
plot(W200, recfft_n(idx));
plot(W200, recfft_sn(idx));
plot(W200, recfft_mn(idx));
plot(W200, mp_fft_n(idx));
plot(W200, mp_sm_fft_n(idx));
plot(W200, FM_fft_n(idx));
plot(W200, FM_mix_fft_n(idx));
plot(W200, spec_fft_n(idx));
legend('True', ...
    sprintf('Kot. (%.2f)',          e1_200), ...
    sprintf('Kot. simps (%.2f)',    e2_200), ...
    sprintf('Mixed (%.2f)',         e3_200), ...
    sprintf('Multipath (%.2f)',     e4_200), ...
    sprintf('MP simps+mix (%.2f)', e5_200), ...
    sprintf('FM (%.2f)',            e6_200), ...
    sprintf('FM mixed (%.2f)',      e7_200), ...
    sprintf('Spectral (%.2f)',      e8_200));
title('FFT Recovery — First 200 Nonneg. Frequencies');
xlabel('\omega'); ylabel('|F(w)| / |F(0)|');

%% Plot 3: Space recovery
f_true = f(grids.D);
f_true_norm = f_true / max(abs(f_true));

rec_a    = align_to_reference(rec,       f_true);
rec_sa   = align_to_reference(rec_simps, f_true);
rec_ma   = align_to_reference(rec_mixed, f_true);
mp_a     = align_to_reference(mp,        f_true);
mp_sma   = align_to_reference(mp_sm,     f_true);
FM_a     = align_to_reference(FM_space,     f_true);
FM_mix_a = align_to_reference(FM_mix_space, f_true);
spec_a   = align_to_reference(spec_space,   f_true);

es1 = rel_linf(rec_a,  f_true_norm');
es2 = rel_linf(rec_sa, f_true_norm');
es3 = rel_linf(rec_ma, f_true_norm');
es4 = rel_linf(mp_a,   f_true_norm');
es5 = rel_linf(mp_sma, f_true_norm');
es6 = rel_linf(FM_a,     f_true_norm');
es7 = rel_linf(FM_mix_a, f_true_norm');
es8 = rel_linf(spec_a,   f_true_norm');

fprintf('\n%-16s  %-12s  %-12s  %-12s\n', 'Algorithm', 'Full FFT', 'FFT (200)', 'Space');
fprintf('%s\n', repmat('-', 1, 56));
fprintf('%-16s  %-12.4f  %-12.4f  %-12.4f\n', 'Kot.',         e1, e1_200, es1);
fprintf('%-16s  %-12.4f  %-12.4f  %-12.4f\n', 'Kot. simps',   e2, e2_200, es2);
fprintf('%-16s  %-12.4f  %-12.4f  %-12.4f\n', 'Mixed',        e3, e3_200, es3);
fprintf('%-16s  %-12.4f  %-12.4f  %-12.4f\n', 'Multipath',    e4, e4_200, es4);
fprintf('%-16s  %-12.4f  %-12.4f  %-12.4f\n', 'MP simps+mix', e5, e5_200, es5);
fprintf('%-16s  %-12.4f  %-12.4f  %-12.4f\n', 'FM',           e6, e6_200, es6);
fprintf('%-16s  %-12.4f  %-12.4f  %-12.4f\n', 'FM mixed',     e7, e7_200, es7);
fprintf('%-16s  %-12.4f  %-12.4f  %-12.4f\n', 'Spectral',     e8, e8_200, es8);
fprintf('%s\n', repmat('-', 1, 56));
%%
figure;
plot(grids.D, f_true_norm,   'k',  'LineWidth', 1.5); hold on;
plot(grids.D, rec_a);
plot(grids.D, rec_sa);
plot(grids.D, rec_ma);
plot(grids.D, mp_a);
plot(grids.D, mp_sma);
plot(grids.D, FM_a);
plot(grids.D, FM_mix_a);
plot(grids.D, spec_a);
legend('True', ...
    sprintf('Kot. (%.2f)',          es1), ...
    sprintf('Kot. simps (%.2f)',    es2), ...
    sprintf('Mixed (%.2f)',         es3), ...
    sprintf('Multipath (%.2f)',     es4), ...
    sprintf('MP simps+mix (%.2f)', es5), ...
    sprintf('FM (%.2f)',            es6), ...
    sprintf('FM mixed (%.2f)',      es7), ...
    sprintf('Spectral (%.2f)',      es8));
title('Space Recovery');
xlabel('x'); ylabel('f(x)');

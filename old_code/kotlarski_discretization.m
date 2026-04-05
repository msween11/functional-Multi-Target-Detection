% Discretization analysis: basic Kotlarski vs Simpson's Kotlarski
% Plots real and imaginary parts of recovered FFT against true FFT

[params, grids] = setparams(n= 10, l = 5, B = 10, z = 1, ...
    t = .01, sigma = 1*10^(0), lambda = .1, r = 1, h = .01);

func_choice = 'f2';
f = choose_function(func_choice);

params.zero_cond = strcmp(func_choice, 'f4');

datastrip = get_obs(params, grids, f);

%% Run both with the same psihat
[ps, psihat, psireg] = get_psihat_auto(params, grids, f, datastrip);
[~, recfft]       = kotlarski(params, grids, psihat, psireg);
[~, recfft_simps] = kotlarski_simps_test(params, grids, psihat, psireg);

%% Reference
mid = length(grids.Wext)/2;
dc  = mid + 1;

true_fft_raw = symfft(f(grids.Dext));

% Normalise by magnitude of DC so both real and imaginary parts are
% on a comparable scale (DC of recovered fft is exp(0)=1 by construction)
norm_c = @(x) x / abs(x(dc));

true_n  = norm_c(true_fft_raw);
rec_n   = norm_c(recfft);
rec_s_n = norm_c(recfft_simps);

%% Errors (L-inf, all normalised by max|true_n| to avoid /~0 when a part is nearly zero)
norm_ref = max(abs(true_n));
comp_err = @(x, ref) max(abs(x - ref)) / norm_ref;

e_re1  = comp_err(real(rec_n),   real(true_n));
e_im1  = comp_err(imag(rec_n),   imag(true_n));
e_all1 = comp_err(rec_n,         true_n);
e_re2  = comp_err(real(rec_s_n), real(true_n));
e_im2  = comp_err(imag(rec_s_n), imag(true_n));
e_all2 = comp_err(rec_s_n,       true_n);

fprintf('\n%-16s  %-12s  %-12s  %-12s\n', 'Algorithm', 'Re L-inf', 'Im L-inf', 'Full L-inf');
fprintf('%s\n', repmat('-', 1, 58));
fprintf('%-16s  %-12.4f  %-12.4f  %-12.4f\n', 'Kotlarski',       e_re1, e_im1, e_all1);
fprintf('%-16s  %-12.4f  %-12.4f  %-12.4f\n', 'Kotlarski simps', e_re2, e_im2, e_all2);
fprintf('%s\n', repmat('-', 1, 58));

%% Plot — full spectrum
figure;

subplot(2,1,1);
plot(grids.Wext, real(true_n),  'k',  'LineWidth', 1.5); hold on;
plot(grids.Wext, real(rec_n),   'b--');
plot(grids.Wext, real(rec_s_n), 'r--');
legend('True', ...
    sprintf('Kotlarski (%.4f)',       e_re1), ...
    sprintf('Kotlarski simps (%.4f)', e_re2));
title('Real part — Full Spectrum');
xlabel('\omega'); ylabel('Re[F(w)]');

subplot(2,1,2);
plot(grids.Wext, imag(true_n),  'k',  'LineWidth', 1.5); hold on;
plot(grids.Wext, imag(rec_n),   'b--');
plot(grids.Wext, imag(rec_s_n), 'r--');
legend('True', ...
    sprintf('Kotlarski (%.4f)',       e_im1), ...
    sprintf('Kotlarski simps (%.4f)', e_im2));
title('Imaginary part — Full Spectrum');
xlabel('\omega'); ylabel('Im[F(w)]');

%% Plot — first 200 nonneg frequencies
idx  = dc : dc + 199;
W200 = grids.Wext(idx);

e_re1_200 = comp_err(real(rec_n(idx)),   real(true_n(idx)));
e_im1_200 = comp_err(imag(rec_n(idx)),   imag(true_n(idx)));
e_re2_200 = comp_err(real(rec_s_n(idx)), real(true_n(idx)));
e_im2_200 = comp_err(imag(rec_s_n(idx)), imag(true_n(idx)));

figure;

subplot(2,1,1);
plot(W200, real(true_n(idx)),  'k',  'LineWidth', 1.5); hold on;
plot(W200, real(rec_n(idx)),   'b--');
plot(W200, real(rec_s_n(idx)), 'r--');
legend('True', ...
    sprintf('Kotlarski (%.4f)',       e_re1_200), ...
    sprintf('Kotlarski simps (%.4f)', e_re2_200));
title('Real part — First 200 Nonneg. Frequencies');
xlabel('\omega'); ylabel('Re[F(w)]');

subplot(2,1,2);
plot(W200, imag(true_n(idx)),  'k',  'LineWidth', 1.5); hold on;
plot(W200, imag(rec_n(idx)),   'b--');
plot(W200, imag(rec_s_n(idx)), 'r--');
legend('True', ...
    sprintf('Kotlarski (%.4f)',       e_im1_200), ...
    sprintf('Kotlarski simps (%.4f)', e_im2_200));
title('Imaginary part — First 200 Nonneg. Frequencies');
xlabel('\omega'); ylabel('Im[F(w)]');

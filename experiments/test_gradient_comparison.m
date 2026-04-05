% test_gradient_comparison.m
%
% Compares two estimates of the horizontal gradient of psihat:
%   fd_gradx    = d14o(psihat, grids.dw, 1)   (4th-order finite differences)
%   psihat_gradx                                (FT of -i*z1*A3, from get_psihat_auto)
%
% The FT-based method is exact by the Fourier derivative property:
%   d/dw1 [ FT[A3](w1,w2) ] = FT[ -i*z1 * A3(z1,z2) ](w1,w2)
% so psihat_gradx is used as the reference. Any discrepancy is finite-
% difference truncation error, which grows with noise level sigma.
%
% Two outputs:
%   Figure 1: diagonal comparison (what Kotlarski actually integrates)
%   Figure 2: relative FD error on diagonal vs noise level sigma

init();
rng(1);

[params, grids] = setparams(n=2^12, l=5, B=10, sigma=0.5, ...
    lambda=0.1, r=0.1, h=0.05);
f = choose_function('f1');
datastrip = get_obs(params, grids, f);

[~, psihat, ~, psihat_gradx] = get_psihat_auto(params, grids, f, datastrip);
fd_gradx = d14o(psihat, grids.dw, 1);

%% --- Diagonal comparison ---
diag_ft = diag(psihat_gradx);
diag_fd = diag(fd_gradx);

rel_err = norm(diag_fd - diag_ft) / norm(diag_ft);
fprintf('sigma = %.2f  |  relative L2 error on diagonal: %.4f\n', params.sigma, rel_err);

% Quick sign check: if the two are strongly anticorrelated, the DFT
% convention uses the opposite sign; flip (-1i*z1) to (+1i*z1) in
% get_psihat_auto and re-run.
r = real(diag_ft)' * real(diag_fd);
fprintf('Correlation of real parts: %.4f  (should be positive)\n', ...
    r / (norm(real(diag_ft)) * norm(real(diag_fd))));

figure('Name', 'Gradient diagonal comparison');
subplot(2,1,1);
plot(grids.Wext, real(diag_ft), 'b', grids.Wext, real(diag_fd), 'r--', 'LineWidth', 1.2);
legend('FT-based', 'Finite diff');
xlabel('\omega');  ylabel('Real part');
title('Re(\partial\psi/\partial\omega_1) on diagonal,  \sigma = 0.5');

subplot(2,1,2);
plot(grids.Wext, imag(diag_ft), 'b', grids.Wext, imag(diag_fd), 'r--', 'LineWidth', 1.2);
legend('FT-based', 'Finite diff');
xlabel('\omega');  ylabel('Imag part');
title('Im(\partial\psi/\partial\omega_1) on diagonal,  \sigma = 0.5');

%% --- Error vs noise level ---
sigmas = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0];
errs   = zeros(size(sigmas));

for s = 1:numel(sigmas)
    [p_s, ~] = setparams(n=2^12, l=5, B=10, sigma=sigmas(s), ...
        lambda=0.1, r=0.1, h=0.05);
    ds_s = get_obs(p_s, grids, f);
    [~, psi_s, ~, gradx_s] = get_psihat_auto(p_s, grids, f, ds_s);
    fd_s  = d14o(psi_s, grids.dw, 1);
    errs(s) = norm(diag(fd_s) - diag(gradx_s)) / norm(diag(gradx_s));
    fprintf('sigma = %.2f  |  relative L2 error: %.4f\n', sigmas(s), errs(s));
end

figure('Name', 'FD error vs noise');
semilogx(sigmas, errs, 'o-k', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
xlabel('\sigma  (noise std)');
ylabel('Relative L2 error  ||fd - FT|| / ||FT||');
title('Finite-difference gradient error relative to FT method');
grid on;

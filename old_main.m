% May 2017
% https://arxiv.org/abs/1705.00641
% https://github.com/NicolasBoumal/MRA

%% Check that Manopt is available
if ~exist('manopt_version', 'file')
    error(sprintf(['Please get Manopt 3.0 or later,\n' ...
                   'available from http://www.manopt.org.\n' ...
                   'Then, run importmanopt.m to add it to the path.'])); %#ok<SPERR>
end

%% Generate a signal x of length N

[params, grids] = setparams(n= 2^10, l = 5, B = 1, z = 1, ...
    t = .01, sigma = .1, lambda = .1, r = .1, h = .05);

func_choice = 'f2';
f = choose_function(func_choice); 

if strcmp(func_choice, 'f4')
   params.zero_cond = 1;
else 
   params.zero_cond = 0;
end

datastrip = get_obs(params, grids, f);

[~, ~, A3M, psihat, psireg] = get_psihat_auto(params, grids, f, datastrip);

b = center_BS(psireg - diag(diag(psireg)));
b2 = center_BS(psireg);
p = 2*grids.L*circshift(max(0,real(diag(psireg))), grids.L/2);

N = length(grids.Dext);
x_true = f(grids.Dext)';

%% Generate M observations with noise variance sigma^2
M = 2^10;
sigma = .1;

tic_generate = tic();
X_data = generate_observations(x_true, M, sigma);
fprintf('Data generation: %.2g [s]\n', toc(tic_generate));

%% Estimate x from the data through invariant features

tic_invfeatmra = tic();

% Estimate the invariants once (this is common to many methods)
tic_invariants = tic();
[mean_est, P_est, B_est] = invariants_from_data(X_data, sigma);
fprintf('Estimation of invariant features: %.2g [s]\n', toc(tic_invariants));

% Estimate the DFT phases from the bispectrum (this can be done with
% different algorithms, and the non-convex approach
% phases_from_bispectrum_real can accept an input z_init to initialize.)
% We could call many different algorithms here to compare.
tic_bispectrum_inversion = tic();
z_est = phases_from_bispectrum_APS_real(b);
fprintf('Estimation of phases from bispectrum: %.2g [s]\n', toc(tic_bispectrum_inversion));

% Recombine (this is computationally cheap)
x_est = combine_features(mean_est, p, z_est);

time_invfeatmra = toc(tic_invfeatmra);

%% Evaluate error, and plot
figure

relerr = relative_error(x_true, x_est); % 2-norm, up to integer shifts
% Align for plotting
x_est = align_to_reference(x_est, x_true);
hold all;
t = 0:(N-1);
plot(t, x_true, '.-');
plot(t, x_est, 'o-');
title(sprintf('Multireference alignment example, M = %d, sigma = %g', M, sigma));
legend('True signal', ...
       sprintf('Invariant features (RMSE: %.2g; time: %.2g [s])', relerr, time_invfeatmra))
      
xlim([0, N-1]);
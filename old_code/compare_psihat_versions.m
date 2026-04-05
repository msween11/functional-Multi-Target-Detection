rng(1);
format compact

func_choice = 'f2';
f = choose_function(func_choice);

[params, grids] = setparams(n=2^10, l=5, B=10, z=1, ...
    t=.001, sigma=.01, lambda=.1, r=.1, h=.05);

params.zero_cond = strcmp(func_choice, 'f4');

datastrip = get_obs(params, grids, f);

% --- compute psihat both ways
tic
% [~, ~, A3M_ref,  psihat_ref,  ~] = get_psihat(params, grids, f, datastrip);
[~, psihat_new,  ~] = get_psihat_auto(params, grids, f, datastrip);
[A3M_ref, psihat_ref,~] = get_psihat_streamed(params, grids, f, datastrip);

toc
%%
% --- errors
A3M_err    = abs(A3M_new   - A3M_ref);
psihat_err = abs(psihat_new - psihat_ref);

fprintf('A3M    rel l2  error: %.4g\n', norm(A3M_err,'fro')   / norm(A3M_ref,'fro'));
fprintf('A3M    rel l-inf error: %.4g\n', max(A3M_err(:))      / max(abs(A3M_ref(:))));
fprintf('psihat rel l2  error: %.4g\n', norm(psihat_err,'fro') / norm(psihat_ref,'fro'));
fprintf('psihat rel l-inf error: %.4g\n', max(psihat_err(:))   / max(abs(psihat_ref(:))));

% --- heatmaps
fig = figure;
tiledlayout(1, 2, 'TileSpacing', 'compact')

nexttile
imagesc(A3M_err)
colorbar
axis square
title('$|A3M_{\rm new} - A3M_{\rm ref}|$', 'Interpreter', 'latex')
xlabel('lag index $j$', 'Interpreter', 'latex')
ylabel('lag index $i$', 'Interpreter', 'latex')

nexttile
imagesc(psihat_err)
colorbar
axis square
title('$|\hat{\psi}_{\rm new} - \hat{\psi}_{\rm ref}|$', 'Interpreter', 'latex')
xlabel('$\omega_1$ index', 'Interpreter', 'latex')
ylabel('$\omega_2$ index', 'Interpreter', 'latex')

dim = [.25 0 .3 .3];
str = sprintf(['$\\sigma = %.3g,\\ \\lambda = %.3g,\\ f = %s,' ...
    '\\ n = 2^{%d}$'], params.sigma, params.lambda, func_choice, log2(params.n));
annotation('textbox', dim, 'String', str, 'Interpreter', 'latex', 'FitBoxToText', 'on');

filename = sprintf('figures/compare_psihat_%s.fig', func_choice);
savefig(fig, filename);

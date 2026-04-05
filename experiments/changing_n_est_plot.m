
format compact
tic

%% ===================== SETUP =====================
seed   = 1:5;
nrange = 5:15;

numSeeds = numel(seed);
numN     = numel(nrange);
numJobs  = numSeeds * numN;

func_choice = 'f2';
f = choose_function(func_choice);

%% ===================== PREALLOC =====================
rel_errs = zeros(numJobs, 1);

dq = parallel.pool.DataQueue;
afterEach(dq, @updateProgress);

streams = parallel.pool.Constant( ...
    @() RandStream('Threefry','Seed',1));

fprintf('Progress: 0 / %d\n', numJobs);

%% ===================== PARFOR =====================
parfor k = 1:numJobs

    [si, li] = ind2sub([numSeeds, numN], k);

    stream = streams.Value;
    stream.Substream = seed(si);
    RandStream.setGlobalStream(stream);

    [params, grids] = setparams(n=2^nrange(li), l=5, B=10, z=1, ...
        t=.001, sigma=.5, lambda=.1, r=.1, h=.05);

    params.zero_cond = strcmp(func_choice, 'f4');

    % --- data strip
    datastrip = get_obs(params, grids, f);

    % --- compute raw moments (no 1/n scaling, needed for n_est formula)
    d = (length(grids.X) - length(grids.x)) / 2;
    datastrip_trunc = datastrip(d:end-d-1);

    data = zeros(length(grids.x), grids.lenD);
    for k2 = 1:grids.lenD
        jj = grids.D(k2) / grids.dx;
        data(:,k2) = datastrip(d+jj:end-d-1+jj);
    end

    A1M = grids.dx * sum(datastrip_trunc);
    A2M = grids.dx * (datastrip_trunc * data);

    % --- n estimator
    covr  = @(x) (params.sigma^2) * exp(-abs(x).^2 ./ (2*params.lambda^2));
    N     = -min(grids.x);
    n_est = A1M^2 / (grids.dx*sum(A2M) - grids.dx*sum(covr(grids.D) .* 2*N));

    rel_errs(k) = abs(n_est - params.n) / params.n;

    send(dq, numJobs);
end

%% ===================== RESHAPE =====================
rel_errs = reshape(rel_errs, [numSeeds, numN]);

%% ===================== MEANS + ERROR BARS =====================
avg_err = mean(rel_errs, 1);
err_std = std(rel_errs, 0, 1);

toc

%% ===================== PLOT =====================
[params, ~] = setparams(n=2^10, l=5, B=10, z=1, ...
    t=.001, sigma=.5, lambda=.1, r=.1, h=.05);

fig = figure;

b = round(polyfit(nrange, log2(avg_err), 1), 2);

errorbar(nrange, log2(avg_err), 1.4427 * err_std ./ avg_err, 'b')

xlabel('$\log_2(n)$', 'Interpreter', 'latex')
ylabel('$\log_2$ Relative Error in $\hat{n}$', 'Interpreter', 'latex')
grid on
axis square
legend(['n\_est  m = ', num2str(b(1))])
title('n estimator convergence vs sample size')

dim = [.25 0 .3 .3];
str = sprintf(['$\\sigma = %.3g,\\ \\lambda = %.3g,\\ f = %s,' ...
    '\\ r = %.3g,\\ h = %.3g$'], ...
    params.sigma, params.lambda, func_choice, params.r, params.h);
annotation('textbox', dim, 'String', str, 'Interpreter', 'latex', 'FitBoxToText', 'on');

filename = sprintf('figures/changing_n_est_%s.fig', func_choice);
savefig(fig, filename);

function updateProgress(N)
    persistent progress
    if isempty(progress)
        progress = 0;
    end
    progress = progress + 1;
    if mod(progress,5)==0 || progress==N
        fprintf('Progress: %d / %d\n', progress, N);
    end
end

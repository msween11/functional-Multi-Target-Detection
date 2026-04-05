format compact
tic

%% ===================== SETUP =====================
seed     = 1:20;
sigrange = [0.001, 0.1:0.1:1];

numSeeds  = numel(seed);
numSigmas = numel(sigrange);
numJobs   = numSeeds * numSigmas;

func_choice  = 'f2';
f = choose_function(func_choice);
show_errorbars = true;

%% ===================== PREALLOC =====================
% 4 algorithms: kott, mpm, fmm, sp
algs = {'kott','mpm','fmm','sp'};
for a = algs
    eval(sprintf('l8_%s_space = zeros(numJobs,1);', a{1}));
    eval(sprintf('l8_%s_freq  = zeros(numJobs,1);', a{1}));
end

dq = parallel.pool.DataQueue;
progress = 0;
afterEach(dq, @updateProgress);

streams = parallel.pool.Constant( ...
    @() RandStream('Threefry','Seed',1));

fprintf('Progress: 0 / %d\n', numJobs);

%% ===================== PARFOR =====================
parfor k = 1:numJobs

    stream = streams.Value;
    stream.Substream = seed(mod(k-1, numSeeds)+1);
    RandStream.setGlobalStream(stream);

    [~, sigi] = ind2sub([numSeeds, numSigmas], k);

    [params, grids] = setparams(n=2^10, l=2, B=10, z=0, ...
        t=.01, sigma=sigrange(sigi), lambda=.1, r=.1, h=.05);
    params.zero_cond = strcmp(func_choice, 'f4');

    datastrip = get_obs(params, grids, f);
    [ps, psihat, psireg, psihat_gradx] = get_psihat_auto(params, grids, f, datastrip);

    % --- all algorithms ---
    [r_kott, f_kott] = kotlarski_auto(params, grids, psihat, psireg, ...
        'method', 'trapz', 'gradx', psihat_gradx);
    [r_mpm,  f_mpm]  = kotlarski_auto(params, grids, psihat, psireg, ...
        'method', 'trapz', 'multipath', true, 'ps', ps, 'gradx', psihat_gradx);
    [r_fmm,  f_fmm]  = freq_march(params, grids, psireg, f, ps);
    [r_sp,   f_sp]   = spectral(params, grids, psireg, f);

    % --- reference ---
    fD  = f(grids.D);
    dc  = length(grids.Wext)/2 + 1;
    tfn = abs(symfft(f(grids.Dext))); tfn = tfn / tfn(dc);

    nf = @(x) abs(x) / abs(x(dc));

    % --- space alignment ---
    r_kott = align_to_reference(r_kott, fD)';
    r_mpm  = align_to_reference(r_mpm,  fD)';
    r_fmm  = align_to_reference(r_fmm,  fD)';
    r_sp   = align_to_reference(r_sp,   fD)';

    % --- freq normalisation ---
    fn_kott = nf(f_kott);
    fn_mpm  = nf(f_mpm);
    fn_fmm  = nf(f_fmm);
    fn_sp   = nf(f_sp);

    % --- errors ---
    rel8s = @(r)  max(abs(r  - fD))  / max(abs(fD));
    rel8f = @(fn) max(abs(fn - tfn)) / max(abs(tfn));

    l8_kott_space(k) = rel8s(r_kott);  l8_kott_freq(k) = rel8f(fn_kott);
    l8_mpm_space(k)  = rel8s(r_mpm);   l8_mpm_freq(k)  = rel8f(fn_mpm);
    l8_fmm_space(k)  = rel8s(r_fmm);   l8_fmm_freq(k)  = rel8f(fn_fmm);
    l8_sp_space(k)   = rel8s(r_sp);    l8_sp_freq(k)   = rel8f(fn_sp);

    send(dq, numJobs);
end
%% ===================== END PARFOR =====================

%% ===================== RESHAPE =====================
r2 = @(x) reshape(x, [numSeeds, numSigmas]);
for a = algs
    eval(sprintf('l8_%s_space = r2(l8_%s_space);', a{1}, a{1}));
    eval(sprintf('l8_%s_freq  = r2(l8_%s_freq);',  a{1}, a{1}));
end

%% ===================== MEANS + STD =====================
avg = @(x) mean(x, 1);
sd  = @(x) std(x, 0, 1) / sqrt(numSeeds);
eb  = @(m,s) 1.4427 * s ./ m;   % delta log2(m) from delta m/m

% pack into structs for cleaner plotting code
for a = algs
    for dom = {'space','freq'}
        vname = sprintf('l8_%s_%s', a{1}, dom{1});
        eval(sprintf('mn.%s = avg(%s);', vname, vname));
        eval(sprintf('er.%s = sd(%s);',  vname, vname));
    end
end

toc

%% ===================== PLOTTING =====================
xax     = sigrange;
xax_log = log2(sigrange);   % log2 kept only for polyfit slope computation

% algorithm display names and colors
names = {'Kot. (T)','MP+mix','FM mix','Spectral'};
cols  = {[0.6 0 0],'g','c','k'};

[params, grids] = setparams(n=2^10, l=5, B=10, z=0, ...
    t=.01, sigma=0.5, lambda=.1, r=.1, h=.05);

dim = [.25 0 .3 .3];
str = sprintf(['$n = 2^{10},\\ \\lambda = %.3g,\\ f = %s,' ...
    '\\ r = %.3g,\\ h = %.3g$'], params.lambda, func_choice, params.r, params.h);

fig = figure;
tiledlayout(1, 2, 'TileSpacing', 'compact');

for tile = 1:2
    nexttile;
    hold on;

    if tile == 1
        dom = 'space';  ttl = 'Changing \sigma, space rel. errors';
    else
        dom = 'freq';   ttl = 'Changing \sigma, freq. rel. errors';
    end

    leg = {};
    for ai = 1:numel(algs)
        a  = algs{ai};
        c  = cols{ai};
        nm = names{ai};

        mn8 = eval(sprintf('mn.l8_%s_%s', a, dom));
        er8 = eval(sprintf('er.l8_%s_%s', a, dom));

        m8 = round(polyfit(xax_log, log2(mn8), 1), 2);

        if show_errorbars
            errorbar(xax, log2(mn8), eb(mn8, er8), 'Color', c, 'LineStyle', '-');
        else
            plot(xax, log2(mn8), 'Color', c, 'LineStyle', '-');
        end

        leg{end+1} = sprintf('%s  m=%.2f', nm, m8(1));
    end

    legend(leg, 'Location', 'northwest', 'FontSize', 6);
    xlabel('$\sigma$', 'Interpreter', 'latex');
    ylabel('$\log_2$ Relative Error', 'Interpreter', 'latex');
    title(ttl);
    grid on;
    axis square;
end

annotation('textbox', dim, 'String', str, 'Interpreter', 'latex', ...
    'FitBoxToText', 'on');

filename = sprintf('figures/changing_sigma_%s.fig', func_choice);
% savefig(fig, filename);

%% ===================== PROGRESS =====================
function updateProgress(N)
    persistent progress
    if isempty(progress), progress = 0; end
    progress = progress + 1;
    if mod(progress, 10) == 0 || progress == N
        fprintf('Progress: %d / %d\n', progress, N);
    end
end

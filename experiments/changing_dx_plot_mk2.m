format compact
tic
%% ===================== SETUP =====================
seed   = 1:4;
lrange = 4:4;

numSeeds = numel(seed);
numL     = numel(lrange);
numJobs  = numSeeds * numL;

func_choice = 'f1';
f = choose_function(func_choice);

%% ===================== PREALLOC (LINEAR) =====================
l8_kot_space  = zeros(numJobs,1);
l8_fm_space   = zeros(numJobs,1);
l8_spec_space = zeros(numJobs,1);

l2_kot_space  = zeros(numJobs,1);
l2_fm_space   = zeros(numJobs,1);
l2_spec_space = zeros(numJobs,1);

l8_kot_freq   = zeros(numJobs,1);
l8_fm_freq    = zeros(numJobs,1);
l8_spec_freq  = zeros(numJobs,1);

l2_kot_freq   = zeros(numJobs,1);
l2_fm_freq    = zeros(numJobs,1);
l2_spec_freq  = zeros(numJobs,1);

dq = parallel.pool.DataQueue;

progress = 0;
afterEach(dq, @updateProgress);

streams = parallel.pool.Constant( ...
    @() RandStream('Threefry','Seed',1));

fprintf('Progress: 0 / %d\n', numJobs);




%% ===================== PARFOR =====================
parfor k = 1:numJobs

    s = seed(k); 
    stream = streams.Value;
    stream.Substream = s;
    RandStream.setGlobalStream(stream);

    % --- flat -> (seed, l)
    [si, li] = ind2sub([numSeeds, numL], k);

    % --- parameters (positional only)
    [params, grids] = setparams(n= 2^10, l = lrange(li),B=10,z=0,...
        t=.01,sigma=.5,lambda = .1,r = .1,h = .05);
    
    params.zero_cond = strcmp(func_choice,'f4');

    % --- data + recon
    datastrip = get_obs(params,grids,f);
    [~,psihat,psireg] = get_psihat_auto(params,grids,f,datastrip);

    [recspace,recfft] = kotlarski_simps_test(params,grids,psihat,psireg);
    [FM_space,FM_freq] = freq_march(params,grids,psireg,f);
    [spec_space,spec_freq] = spectral(params,grids,psireg,f);

    % --- reference (evaluate once)
    fD = f(grids.D);

    % --- alignment
    recspace = align_to_reference(recspace,fD)';

    plotting = real(symfft(f(grids.Dext)));
    plotting = plotting ./ plotting(numel(plotting)/2 + 1);

    fmn = FM_freq   ./ FM_freq(numel(FM_freq)/2 + 1);
    sn  = spec_freq ./ spec_freq(numel(spec_freq)/2 + 1);

    % --- errors (ONLY index k)
    l8_kot_space(k)  = max(abs(recspace - fD)) / max(abs(fD));
    l8_fm_space(k)   = max(abs(FM_space - fD)) / max(abs(fD));
    l8_spec_space(k) = max(abs(spec_space - fD)) / max(abs(fD));

    l2_kot_space(k)  = norm(recspace - fD) / norm(fD);
    l2_fm_space(k)   = norm(FM_space - fD) / norm(fD);
    l2_spec_space(k) = norm(spec_space - fD) / norm(fD);

    l8_kot_freq(k)   = max(abs(real(recfft) - plotting)) / max(abs(plotting));
    l8_fm_freq(k)    = max(abs(real(fmn)   - plotting)) / max(abs(plotting));
    l8_spec_freq(k)  = max(abs(real(sn)    - plotting)) / max(abs(plotting));

    l2_kot_freq(k)   = norm(real(recfft) - plotting) / norm(plotting);
    l2_fm_freq(k)    = norm(real(fmn)   - plotting) / norm(plotting);
    l2_spec_freq(k)  = norm(real(sn)    - plotting) / norm(plotting);

    send(dq, numJobs);
end
%% ===================== END PARFOR =====================

%% ===================== RESHAPE =====================
reshape2 = @(x) reshape(x,[numSeeds,numL]);

l8_kot_space  = reshape2(l8_kot_space);
l8_fm_space   = reshape2(l8_fm_space);
l8_spec_space = reshape2(l8_spec_space);

l2_kot_space  = reshape2(l2_kot_space);
l2_fm_space   = reshape2(l2_fm_space);
l2_spec_space = reshape2(l2_spec_space);

l8_kot_freq   = reshape2(l8_kot_freq);
l8_fm_freq    = reshape2(l8_fm_freq);
l8_spec_freq  = reshape2(l8_spec_freq);

l2_kot_freq   = reshape2(l2_kot_freq);
l2_fm_freq    = reshape2(l2_fm_freq);
l2_spec_freq  = reshape2(l2_spec_freq);

%% ===================== MEANS + ERROR BARS =====================
avg = @(x) mean(x,1);
std_dev = @(x) std(x,0,1);

avg_table = table( ...
    avg(l8_kot_space)',  avg(l8_fm_space)',  avg(l8_spec_space)', ...
    avg(l2_kot_space)',  avg(l2_fm_space)',  avg(l2_spec_space)', ...
    avg(l8_kot_freq)',   avg(l8_fm_freq)',   avg(l8_spec_freq)', ...
    avg(l2_kot_freq)',   avg(l2_fm_freq)',   avg(l2_spec_freq)', ...
    'VariableNames', { ...
    'l8_kot_space_err','l8_fm_space_err','l8_spec_space_err', ...
    'l2_kot_space_err','l2_fm_space_err','l2_spec_space_err', ...
    'l8_kot_freq_err','l8_fm_freq_err','l8_spec_freq_err', ...
    'l2_kot_freq_err','l2_fm_freq_err','l2_spec_freq_err'});
errs = table( ...
    std_dev(l8_kot_space)',  std_dev(l8_fm_space)',  std_dev(l8_spec_space)', ...
    std_dev(l2_kot_space)',  std_dev(l2_fm_space)',  std_dev(l2_spec_space)', ...
    std_dev(l8_kot_freq)',   std_dev(l8_fm_freq)',   std_dev(l8_spec_freq)', ...
    std_dev(l2_kot_freq)',   std_dev(l2_fm_freq)',   std_dev(l2_spec_freq)', ...
    'VariableNames', avg_table.Properties.VariableNames);

toc

[params, grids] = setparams(n= 2^10, l = 5,B=10,z=0,...
   t=.01,sigma=.1,lambda = .1,r = .1,h = .05);
    
%% now we will make our plot(s) and export the .fig file


fig = figure
tiledlayout(1,2, 'TileSpacing', 'compact')
%tile 1
nexttile

%slopes of best fit lines
b8_kot = round(polyfit(lrange, log2(avg_table.l8_kot_space_err),1),2);
b8_fm = round(polyfit(lrange, log2(avg_table.l8_fm_space_err),1),2);
b8_spec = round(polyfit(lrange, log2(avg_table.l8_spec_space_err),1),2);

b2_kot = round(polyfit(lrange, log2(avg_table.l2_kot_space_err),1),2);
b2_fm = round(polyfit(lrange, log2(avg_table.l2_fm_space_err),1),2);
b2_spec = round(polyfit(lrange, log2(avg_table.l2_spec_space_err),1),2);

%recording parameters
dim = [.25 0 .3 .3];
str = sprintf(['$\\sigma = %.3g,\\ \\lambda = %.3g,\\ f = %s,' ...
    '\\ r = %.3g,\\ h = %.3g$'],...
    params.sigma, params.lambda, func_choice, params.r, params.h);
annotation('textbox',dim,'String',str, 'Interpreter', 'latex',...
    'FitBoxToText','on');

%plotting l8 errs
errorbar(lrange, log2(avg_table.l8_kot_space_err), ...
    1.4427*errs.l8_kot_space_err ./ avg_table.l8_kot_space_err, 'red')
hold on
errorbar(lrange, log2(avg_table.l8_fm_space_err), ...
    1.4427*errs.l8_fm_space_err ./ avg_table.l8_fm_space_err,'blue')
errorbar(lrange, log2(avg_table.l8_spec_space_err), ...
    1.4427*errs.l8_spec_space_err ./ avg_table.l8_spec_space_err,'green')

%plotting l2 errs
errorbar(lrange, log2(avg_table.l2_kot_space_err), ...
    1.4427*errs.l2_kot_space_err ./ avg_table.l2_kot_space_err, 'Color','red','LineStyle' ,'--')
hold on
errorbar(lrange, log2(avg_table.l2_fm_space_err), ...
    1.4427*errs.l2_fm_space_err ./ avg_table.l2_fm_space_err,'Color','blue','LineStyle' ,'--')
errorbar(lrange, log2(avg_table.l2_spec_space_err), ...
    1.4427*errs.l2_spec_space_err ./ avg_table.l2_spec_space_err,'Color','green','LineStyle' ,'--')

%aesthetics
xlabel('$-\log_2(dx)$','interpreter','latex')
ylabel([' $$\log_2 $$ Relative $$l^\infty $$ and $$l^2 $$ Errors'],...
    'Interpreter','latex')
grid on
axis square
legend(['Kotlarski l8' ' m = ', num2str(b8_kot(1))],['FM l8' ' m = ', num2str(b8_fm(1))],...
    ['Spectral l8' ' m = ', num2str(b8_spec(1))], ['Kotlarski l2' ' m = ', num2str(b2_kot(1))],...
    ['FM l2' ' m = ', num2str(b2_fm(1))],...
    ['Spectral l2' ' m = ', num2str(b2_spec(1))])
title('Changing dx, space rel errors')

%tile 2----------------------------------------------
nexttile

%slopes of best fit lines
b8_kot = round(polyfit(lrange, log2(avg_table.l8_kot_freq_err),1),2);
b8_fm = round(polyfit(lrange, log2(avg_table.l8_fm_freq_err),1),2);
b8_spec = round(polyfit(lrange, log2(avg_table.l8_spec_freq_err),1),2);

b2_kot = round(polyfit(lrange, log2(avg_table.l2_kot_freq_err),1),2);
b2_fm = round(polyfit(lrange, log2(avg_table.l2_fm_freq_err),1),2);
b2_spec = round(polyfit(lrange, log2(avg_table.l2_spec_freq_err),1),2);

%plotting l8 errs
errorbar(lrange, log2(avg_table.l8_kot_freq_err), ...
    1.4427*errs.l8_kot_freq_err ./ avg_table.l8_kot_freq_err, 'red')
hold on
errorbar(lrange, log2(avg_table.l8_fm_freq_err), ...
    1.4427*errs.l8_fm_freq_err ./ avg_table.l8_fm_freq_err,'blue')
errorbar(lrange, log2(avg_table.l8_spec_freq_err), ...
    1.4427*errs.l8_spec_freq_err ./ avg_table.l8_spec_freq_err,'green')

%plotting l2 errs
errorbar(lrange, log2(avg_table.l2_kot_freq_err), ...
    1.4427*errs.l2_kot_freq_err ./ avg_table.l2_kot_freq_err, 'Color','red','LineStyle' ,'--')
hold on
errorbar(lrange, log2(avg_table.l2_fm_freq_err), ...
    1.4427*errs.l2_fm_freq_err ./ avg_table.l2_fm_freq_err,'Color','blue','LineStyle' ,'--')
errorbar(lrange, log2(avg_table.l2_spec_freq_err), ...
    1.4427*errs.l2_spec_freq_err ./ avg_table.l2_spec_freq_err,'Color','green','LineStyle' ,'--')

%aesthetics
xlabel('$-\log_2(dx)$','interpreter','latex')
ylabel([' $$\log_2 $$ Relative $$l^\infty $$ and $$l^2 $$ Errors'],...
    'Interpreter','latex')
grid on
axis square
legend(['Kotlarski l8' ' m = ', num2str(b8_kot(1))],['FM l8' ' m = ', num2str(b8_fm(1))],...
    ['Spectral l8' ' m = ', num2str(b8_spec(1))], ['Kotlarski l2' ' m = ', num2str(b2_kot(1))],...
    ['FM l2' ' m = ', num2str(b2_fm(1))],...
    ['Spectral l2' ' m = ', num2str(b2_spec(1))])
title('Changing dx, freq. rel errors')

%saving the figure
filename = sprintf('figures/changing_dx_%s.fig', func_choice); 
savefig(fig, filename);

function updateProgress(N)
    persistent progress
    if isempty(progress)
        progress = 0;
    end
    progress = progress + 1;

    % throttle output (important)
    if mod(progress,10)==0 || progress==N
        fprintf('Progress: %d / %d\n', progress, N);
    end
end
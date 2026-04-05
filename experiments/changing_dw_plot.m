
%to do: ADD l2 error, freq errors,
% ADD: manifold opt + APS init by both kot and FM

format compact
seed = 1:20;
Brange = 1:20;
numSeeds = numel(seed);
tables = cell(numSeeds,1);

func_choice = 'f2';
f = choose_function(func_choice); 



streams = parallel.pool.Constant( ...
    @() RandStream('Threefry','Seed',1));

parfor si = 1:numSeeds

    %set seed
    s = seed(si); 
    stream = streams.Value;
    stream.Substream = s;
    RandStream.setGlobalStream(stream);

    %clear and set data vectors
    
    %l8 space errs
    l8_kot_space_err = zeros(1,length(Brange));
    l8_fm_space_err = zeros(1,length(Brange));
    l8_spec_space_err = zeros(1,length(Brange));

    %l2 space errs
    l2_kot_space_err = zeros(1,length(Brange));
    l2_fm_space_err = zeros(1,length(Brange));
    l2_spec_space_err = zeros(1,length(Brange));

    %l8 freq errs
    l8_kot_freq_err = zeros(1,length(Brange));
    l8_fm_freq_err = zeros(1,length(Brange));
    l8_spec_freq_err = zeros(1,length(Brange));

    %l2 freq errs
    l2_kot_freq_err = zeros(1,length(Brange));
    l2_fm_freq_err = zeros(1,length(Brange));
    l2_spec_freq_err = zeros(1,length(Brange));



    %start Brange loop
    for idx = 1:numel(Brange)
    
    %initialize parameters and grids
    [params, grids] = setparams(n = 2^10, l = 5, B = Brange(idx), z = 1, ...
    t = .001, sigma = .1, lambda = .1, r = .1, h = .05);
    
    %double check zero_cond for flipping or not
    if strcmp(func_choice, 'f4')
        params.zero_cond = 1;
    else 
        params.zero_cond = 0;
    end
   
    %actual algorithm code
    datastrip = get_obs(params, grids, f);
    [~, psihat,psireg] = get_psihat_auto(params, grids, f, datastrip);
    [recspace, recfft] = kotlarski_simps_test(params, grids, psihat, psireg);
    [FM_space, FM_freq] = freq_march(params, grids, psireg,f);
    [spectral_space, spectral_freq] =  spectral(params,grids, psireg, f);


    %aligning and preparing to collect errors
    recspace = align_to_reference(recspace, f(grids.D))';
    plotting = real(symfft(f(grids.Dext)));
    plotting = plotting ./plotting(length(grids.Dext)/2+1);
    fmn = FM_freq ./ FM_freq(length(grids.Dext)/2+1);
    sn = spectral_freq ./ spectral_freq(length(grids.Dext)/2+1);

    %collect errors
    l8_kot_space_err(idx) = max(abs(recspace-f(grids.D)))/max(abs(f(grids.D)));
    l8_fm_space_err(idx) = max(abs(FM_space-f(grids.D)))/max(abs(f(grids.D)));
    l8_spec_space_err(idx) = max(abs(spectral_space-f(grids.D)))/max(abs(f(grids.D)));

    l2_kot_space_err(idx) = norm(recspace-f(grids.D))/norm(f(grids.D));
    l2_fm_space_err(idx) = norm(FM_space-f(grids.D))/norm(f(grids.D));
    l2_spec_space_err(idx) = norm(spectral_space-f(grids.D))/norm(f(grids.D));

    l8_kot_freq_err(idx) = max(abs(real(recfft)-plotting))/max(abs(plotting));
    l8_fm_freq_err(idx) = max(abs(real(fmn)-plotting))/max(abs(plotting));
    l8_spec_freq_err(idx) = max(abs(real(sn)-plotting))/max(abs(plotting));

    l2_kot_freq_err(idx) = norm(real(recfft)-plotting)/norm(plotting);
    l2_fm_freq_err(idx) = norm(real(fmn)-plotting)/norm(plotting);
    l2_spec_freq_err(idx) = norm(real(sn)-plotting)/norm(plotting);

    end

%record data for s seed table
T = table(l8_kot_space_err', l8_fm_space_err', l8_spec_space_err', ...
    l2_kot_space_err', l2_fm_space_err', l2_spec_space_err',...
    l8_kot_freq_err', l8_fm_freq_err', l8_spec_freq_err',...
    l2_kot_freq_err', l2_fm_freq_err', l2_spec_freq_err');
T.Properties.VariableNames =["l8_kot_space_err", "l8_fm_space_err", "l8_spec_space_err", ...
    "l2_kot_space_err", "l2_fm_space_err", "l2_spec_space_err",...
    "l8_kot_freq_err", "l8_fm_freq_err", "l8_spec_freq_err",...
    "l2_kot_freq_err", "l2_fm_freq_err", "l2_spec_freq_err"];
tables{si} = T;

%end of rng loop below

end

%after the rng loops, collect the final data table:
avg_table = tables{1};
for si = 2:numSeeds
    avg_table = avg_table + tables{si};
end
avg_table = avg_table ./ numSeeds;


%now we will calculate the error bars
errs = 0;
for si = 1:numSeeds
   errs = errs + (tables{si} - avg_table).^2;
end

errs = 1/(length(seed)-1).*errs;
errs = sqrt(errs);
%errs is the final table for error bars

[params, grids] = setparams(n = 2^10, l = 5, B = 10, z = 1, ...
    t = .001, sigma = .1, lambda = .1, r = .1, h = .05);
    

%% now we will make our plot(s) and export the .fig file


fig = figure
tiledlayout(1,2, 'TileSpacing', 'compact')
%tile 1
nexttile

%slopes of best fit lines
b8_kot = round(polyfit(Brange, log2(avg_table.l8_kot_space_err),1),2);
b8_fm = round(polyfit(Brange, log2(avg_table.l8_fm_space_err),1),2);
b8_spec = round(polyfit(Brange, log2(avg_table.l8_spec_space_err),1),2);

b2_kot = round(polyfit(Brange, log2(avg_table.l2_kot_space_err),1),2);
b2_fm = round(polyfit(Brange, log2(avg_table.l2_fm_space_err),1),2);
b2_spec = round(polyfit(Brange, log2(avg_table.l2_spec_space_err),1),2);

%recording parameters
dim = [.25 0 .3 .3];
str = sprintf(['$\\sigma = %.3g,\\ \\lambda = %.3g,\\ f = %s,' ...
    '\\ r = %.3g,\\ h = %.3g$'],...
    params.sigma, params.lambda, func_choice, params.r, params.h);
annotation('textbox',dim,'String',str, 'Interpreter', 'latex',...
    'FitBoxToText','on');

%plotting l8 errs
errorbar(Brange, log2(avg_table.l8_kot_space_err), ...
    1.4427*errs.l8_kot_space_err ./ avg_table.l8_kot_space_err, 'red')
hold on
errorbar(Brange, log2(avg_table.l8_fm_space_err), ...
    1.4427*errs.l8_fm_space_err ./ avg_table.l8_fm_space_err,'blue')
errorbar(Brange, log2(avg_table.l8_spec_space_err), ...
    1.4427*errs.l8_spec_space_err ./ avg_table.l8_spec_space_err,'green')

%plotting l2 errs
errorbar(Brange, log2(avg_table.l2_kot_space_err), ...
    1.4427*errs.l2_kot_space_err ./ avg_table.l2_kot_space_err, 'Color','red','LineStyle' ,'--')
hold on
errorbar(Brange, log2(avg_table.l2_fm_space_err), ...
    1.4427*errs.l2_fm_space_err ./ avg_table.l2_fm_space_err,'Color','blue','LineStyle' ,'--')
errorbar(Brange, log2(avg_table.l2_spec_space_err), ...
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
b8_kot = round(polyfit(Brange, log2(avg_table.l8_kot_freq_err),1),2);
b8_fm = round(polyfit(Brange, log2(avg_table.l8_fm_freq_err),1),2);
b8_spec = round(polyfit(Brange, log2(avg_table.l8_spec_freq_err),1),2);

b2_kot = round(polyfit(Brange, log2(avg_table.l2_kot_freq_err),1),2);
b2_fm = round(polyfit(Brange, log2(avg_table.l2_fm_freq_err),1),2);
b2_spec = round(polyfit(Brange, log2(avg_table.l2_spec_freq_err),1),2);

%plotting l8 errs
errorbar(Brange, log2(avg_table.l8_kot_freq_err), ...
    1.4427*errs.l8_kot_freq_err ./ avg_table.l8_kot_freq_err, 'red')
hold on
errorbar(Brange, log2(avg_table.l8_fm_freq_err), ...
    1.4427*errs.l8_fm_freq_err ./ avg_table.l8_fm_freq_err,'blue')
errorbar(Brange, log2(avg_table.l8_spec_freq_err), ...
    1.4427*errs.l8_spec_freq_err ./ avg_table.l8_spec_freq_err,'green')

%plotting l2 errs
errorbar(Brange, log2(avg_table.l2_kot_freq_err), ...
    1.4427*errs.l2_kot_freq_err ./ avg_table.l2_kot_freq_err, 'Color','red','LineStyle' ,'--')
hold on
errorbar(Brange, log2(avg_table.l2_fm_freq_err), ...
    1.4427*errs.l2_fm_freq_err ./ avg_table.l2_fm_freq_err,'Color','blue','LineStyle' ,'--')
errorbar(Brange, log2(avg_table.l2_spec_freq_err), ...
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
filename = sprintf('figures/changing_dw_%s.fig', func_choice); 
savefig(fig, filename);
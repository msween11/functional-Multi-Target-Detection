function datastrip = get_obs(params, grids,f, options)

arguments
    params struct
    grids struct
    f function_handle
    options.type = '';
end

%creating local shifts
shifts = zeros(1, params.n);
shifts(1) = 0;
local_shifts = params.sd*rand(1,params.n);

%creating global shifts, well separated
for j = params.n/2+1:params.n
    shifts(j) = shifts(j-1)+4+grids.dx+local_shifts(j);
end

for j = params.n/2:-1:1
    shifts(j) = shifts(j+1)-4-grids.dx-local_shifts(j);
end

if strcmp(options.type, 'discrete')
    %DISCRETE SETTING

    shifts_on_grid = grids.dx * floor(shifts / grids.dx);

    %reindexing to start from 1
    centers = round((shifts_on_grid - grids.X(1)) / grids.dx) + 1;
    func = f(grids.D); dd = (length(func))/2;

    M = zeros(1,length(grids.X));
    for j = 1:params.n
        c = centers(j);
        M(c-dd:c+dd-1) = M(c-dd:c+dd-1) + func;
        % M = M + circshift(f(grids.X), )
    end

    
else 
    %FUNCTIONAL SETTING

    %creating M, noiseless observation
    M = zeros(1,length(grids.X));
    for j = 2:params.n
        start = floor((shifts(j)-1) / grids.dx) * grids.dx;  
        stop = ceil((shifts(j)+1) / grids.dx) * grids.dx;  
        idx_vals = start:grids.dx:stop;
        idx = round((idx_vals - grids.X(1)) / grids.dx) + 1;
        M(idx) = M(idx) + f(idx_vals - shifts(j));
    end
end

%adding noise in now to create datastrip
ep = fast_gp_nonperiodic(length(grids.X), grids.dx, params.sigma, params.lambda);
datastrip = M + ep;


end
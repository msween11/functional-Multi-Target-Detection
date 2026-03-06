function [params, grids] = setparams(opts)

arguments
    opts.n (1,1) double = 1e5
    opts.l (1,1) double = 5
    opts.B (1,1) double = 10
    opts.z (1,1) double = 0
    opts.t (1,1) double = .001
    opts.sigma (1,1) double = 1
    opts.lambda (1,1) double = .1
    opts.r (1,1) double = .01
    opts.h (1,1) double = .05
    opts.ls (1,1) double = 0.5;
end

params = struct();
params.n = opts.n; 
params.l = opts.l; 
params.B = opts.B;
params.zero_cond = opts.z; 
params.thresh = opts.t; 
params.sigma = opts.sigma; 
params.lambda = opts.lambda;
params.r = opts.r; 
params.h = opts.h;
params.sd = .5;

N = ceil((2+params.sd)*(params.n+1)); 
dx = 2^(-params.l);

grids = struct();
grids.dx = 2^(-params.l);
grids.x = -N:dx:N-dx;
grids.X = -N-3:dx:N+3-dx;
grids.D = -2:dx:2-dx;
grids.Dext = -3*params.B:dx:3*params.B-dx;
grids.Wext = -pi/dx:(pi/(3*params.B)):pi/dx-(pi/(3*params.B)); 
grids.lenD = length(grids.D);
grids.dw = grids.Wext(2) - grids.Wext(1);
grids.L = length(grids.Wext);

end
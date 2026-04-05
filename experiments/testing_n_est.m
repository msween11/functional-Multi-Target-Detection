rng(1);
[params, grids] = setparams(n= 2^10, l = 5, B = 10, z = 1, ...
    t = .001, sigma = .5,lambda = .1, r = .1, h = .05);

func_choice = 'f2';
f = choose_function(func_choice); 

if strcmp(func_choice, 'f4')
   params.zero_cond = 1;
else 
   params.zero_cond = 0;
end

sd = .5;
shift_centers = zeros(1, params.n);
shift_centers(1) = 0;

for j = params.n/2+1:params.n
    shift_centers(j) = shift_centers(j-1)+4+sd*rand;
end

for j = params.n/2:-1:1
    shift_centers(j) = shift_centers(j+1)-4-sd*rand;
end

M = zeros(1,length(grids.X));
func = f(-2:grids.dx:2);
dd = (length(func)-1)/2;

centers = round((shift_centers - grids.X(1)) / grids.dx) + 1;

for j = 1:params.n 
    c = centers(j);
    M(c-dd:c+dd) = M(c-dd:c+dd) + func;
end

ep = fast_gp_nonperiodic(length(grids.X), grids.dx, params.sigma, params.lambda);
datastrip = M + ep;

d = length(grids.X)-length(grids.x); d=d/2;
datastrip_trunc = datastrip(d:end-d-1);

    data = zeros(length(grids.x), length(grids.D)); %data for A3M
    
    for k = 1:grids.lenD
        j = grids.D(k)/grids.dx;
        data(:,k) = datastrip(d+j:end-d-1+j);
    end

    A1M =  grids.dx * sum(datastrip_trunc);
    A2M =  grids.dx * (datastrip_trunc * data);   

covr = @(x) (params.sigma^2)*exp(-abs(x).^2./(2*params.lambda^2)); 
u = covr(grids.D)*integral(f,-1,1);
%%
epdata = zeros(length(grids.x), length(grids.D));
ep_trunc = ep(d:end-d-1);

for k = 1:grids.lenD
    j = grids.D(k)/grids.dx;
    epdata(:,k) = ep(d+j:end-d-1+j);
end
N = -min(grids.x);
A2ep =  grids.dx * (ep_trunc * epdata);  

n_est = A1M^2./((grids.dx*sum(A2M))-grids.dx*sum(covr(grids.D).*2*N));

err = norm(n_est-params.n)./norm(params.n)



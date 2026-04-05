tic
rng(1);
n = 2^10; 
N = 3*(n+1); l = 5; B = 10;
dx = 2^(-l);
x=-N:2^(-l):N-2^(-l);
X=-N-3:2^(-l):N+3-2^(-l);
D = -2:2^(-l):2-2^(-l);
Dext = -3*B:2^(-l):3*B-2^(-l);
Wext = -pi*(2^l):(pi/(3*B)):pi*(2^l)-(pi/(3*B)); 
dw = Wext(2) - Wext(1);
L = length(Wext);
lenD = length(D);
zero_cond = 0;
thresh = 0.001;
sigma = .000001; lambda = .1;
r = 1; h = 0.041;

[params, grids ] = setparams(  10^5,5,10,0,.001, .000001,.1,.01,.05);
% ...check Anna's approx of n

symfft = @(x) ifftshift(fft(fftshift(x)));

% LINEAR SPLINE
% f1 = @(x) max(1-abs(x),0);
% f1 = @(x) f1(x) + f1(pi*x); 
% nrm = @(x) max(abs(f1(x)));
% f = @(x) f1(x) / nrm(x);

% CUBIC SPLINE
f2 = @(x) (1/6).*( 0.* (x<0)+...
    (x.^3) .* (0<= x & x < 1)+...
    (-3*(x-1).^3 + 3*(x-1).^2+3*(x-1)+1) .* (1<= x & x < 2)+...
    (3*(x-2).^3-6*(x-2).^2+4) .* (2 <= x & x < 3)+...
    ((4-x).^3) .* (3 <=x & x < 4)+...
    0 .* (x >= 4));
f2 = @(x) f2(2.*x+2);
f2 = @(x) f2(x) + f2(pi.*x);
nrm = @(x) max(abs(f2(x)));
f = @(x) f2(x) ./ nrm(x);

% SHIFTED GABOR
% k=8;
% f3 = @(x) exp(-20*x.^2).*cos(k*x);
% f3 = @(x) f3(x-.3);
% nrm = @(x) max(abs(f3(x)));
% f = @(x) f3(x) / nrm(x);

% GMM2, shifted
% means = [-.4 ;.4]; devs = [.02];
% f4 = gmdistribution(means, devs);
% f4 = @(x) pdf(f4, x')';
% nrm = @(x) max(abs(f4(x)));
% f = @(x) f4(x) / nrm(x);
% f = @(x) f(x-.1);


% SAMPLING SHIFTS
minsep = 2;
s = -minsep + 2*minsep*rand(1, n); s(1) = abs(s(1)); 
shifts = zeros(1,n);

for k = 2:2:n
    shifts(k) = -3*k + s(k-1);
    if k >2 && abs(shifts(k) - shifts(k-2)) < minsep
        shifts(k) = shifts(k-2) - minsep * sign(shifts(k) - shifts(k-2));
    end
end

for k = 3:2:n
    shifts(k) = 3*k + s(k-1);
    if abs(shifts(k) - shifts(k-2)) < minsep    
        shifts(k) = shifts(k-2) + minsep * sign(shifts(k) - shifts(k-2));
    end
end

%  CONSTRUCTING M(t)
M = zeros(1,length(X));
origin = length(X)/2;
func = f(-2:2^(-l):2);
dd = (length(func)-1)/2;
for j = 1:n 
    c = origin + round(shifts(j)*(2^l));
    M(c-dd:c+dd) = M(c-dd:c+dd) + func;
end

% CONSTRUCTING NOISE 
ep = fast_gp_nonperiodic(length(X), 2^(-l), sigma, lambda);
datastrip = M + ep;


%% constructing A3M
tic
d = length(X)-length(x); d=d/2;
dd = lenD/2;
datastrip_trunc = datastrip(d:end-d-1);
try
    data = zeros(length(x), length(D)); %data for A3M
    
    for k = 1:length(D)
        j = D(k)*2^l;
        data(:,k) = datastrip(d+j:end-d-1+j);
    end

    A1M =  2^(-l) * sum(datastrip_trunc);
    A2M =  2^(-l) * (datastrip_trunc * data);   
    A3M = (1/n)*2^(-l) * ((data .* datastrip_trunc')' * data); 
catch ME
    if strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded')
        fprintf('Vectorized method exceeded memory limits — switching to streaming.\n');
        %% fully streamed approach.  
                
        UL = zeros(dd+1,dd+1); UR = zeros(dd+1,dd-1); 
        DD = D(1:dd+1); DD2 = D(dd+2:end);                
                
        for i = 1:dd+1
            ii = DD(i)*(2^l);
            di = datastrip_trunc .* datastrip(d+ii:end-d-1+ii);
            for  j = 1:i
                jj = DD(j)*2^l;
                UL(i,j) = 1/n*2^(-l)*dot(di,datastrip(d+jj:end-d-1+jj));
                UL(j,i) = UL(i,j);
            end

            for  j = 1:dd-1
                jj = DD2(j)*(2^l);
                UR(i,j) = 1/n*2^(-l)*dot(di,datastrip(d+jj:end-d-1+jj));
            end
        end
        LR = flip(flip(UL,1),2);%rot90(UL,2);
        LR = LR(2:end-1, 2:end-1);

        A3Ms = [UL UR ; UR' LR]  ;   
    end
end
fprintf('A3M constructed. \n');
toc

%%  UNBIASING  

covr = @(x) (sigma^2)*exp(-abs(x).^2./(2*lambda^2)); 
u = covr(D)*integral(f,-1,1);
A3M = A3M-u-u'-(integral(f,-1,1)*covr(D-D'));

%%

n_est = A1M^2./((dx*sum(A2M))-4*dx*sum(covr(D)))
%% PADDING A3M AND MAPPING INTO FREQUENCY 

n1 = length(D);
n2 = length(Dext);
nside = (n2-n1)/2;
A3M_padded = [zeros(nside,n2); zeros(n1,nside) A3M zeros(n1,nside); zeros(nside,n2)];

F = dftmtx(length(Dext));
A3M_ft = center_BS(F'*center_BS(A3M_padded)*F);


%% DOING KOTLARSKI NOW

psihat = (1/L)*A3M_ft;
% Doing Kotlarski's identity EMPIRICALLY

[empx, empy] = gradient(psihat, dw);

adx = diag(empx); 
ad_emp = diag(psihat);

integrand_emp = adx./ad_emp; integrand_emp = integrand_emp';

rec_emp = zeros(1,length(Wext)/2);
rec_emp(1) = 0;
for j = 2:length(Wext)/2
    rec_emp(j) = rec_emp(j-1) + dw*integrand_emp( length(Wext)/2 +j);
end
rec_emp = exp(rec_emp);

psireg = psihat./(L*min(ones(length(Wext))/L, r*sqrt(n)*abs(psihat)));

ad_reg = diag(psireg);
integrand_reg = adx./ad_reg; integrand_reg = integrand_reg';

rec_reg = zeros(1,length(Wext)/2);
rec_reg(1) = 0;
for j = 2:length(Wext)/2
    rec_reg(j) = rec_reg(j-1) + dw*integrand_reg( length(Wext)/2 +j);
end

rec1_reg = (exp(rec_reg)); 

if(zero_cond == 1)
    growth_filt = .001;
    T = thresh;
    ps = real(ad_emp)';
    g = gradient(ps, dw);
    gg = gradient(g, dw);

    pidz = findz(g);
    pidz_filt = pidz(abs(gg(pidz))> growth_filt );
    pidz_filt = pidz_filt(Wext(pidz_filt)>0);
    pz_filt = Wext(pidz_filt);
    fliploc = pz_filt(abs(ps(pidz_filt))< T); %this is by thresh

    %getting the index
    [~,loc] = ismember(fliploc, Wext);
    loc = loc-length(Wext)/2;
    loc = [loc length(Wext)/2];

    %flipping
    for j = 1:length(loc)-1
        rec1_reg(loc(j)+1:loc(j+1))=(-1)^j*rec1_reg(loc(j)+1:loc(j+1));
    end

    %zeroing the little windows around the zeros
    ep = 2; zloc = [];
    for w=1:ep
        temp = [loc-w loc loc+w];
        zloc = [zloc temp];
    end
    clear('temp');
    zloc = unique(zloc);
    zloc(end-ep+1:end) = [];     
    rec1_reg(zloc) = 0;
end

plotting = symfft(f(Dext));
plotting = plotting ./plotting(length(Dext)/2+1);
recfft = [conj(flip(rec1_reg)) (rec1_reg)];
l8_freq_rel_err = max(abs(recfft - plotting))/max(abs(plotting));

%hfix recovery in space now
phiK = @(w) 1 .* (-1 <= w & w <= 1);


rec = zeros(1, length(Dext));
for j = 1:length(Dext)
    rec(j) = (dw/(2*pi)).* sum(exp(1i.*Wext*Dext(j)).*...
        recfft .* phiK(h*Wext));
end
t = length(D)/2;  
tt = length(Dext)/2;
rec_trunc = real(rec(tt-t+1:tt+t)); %choosing real here should be ok
recspace = rec_trunc/max(abs(rec_trunc));

%

Dd = length(Dext)/2;
p = circshift(real(diag(psireg)), Dd); %choose psihat here maybe?

phases = phases_from_bispectrum_FM_real(center_BS(psireg));
FM_freq = sqrt(p).*phases;
FM_freq = circshift(FM_freq', Dd);

FM_space = zeros(1, length(Dext));
for j = 1:length(Dext)
    FM_space(j) = (dw/(2*pi)).* sum(exp(1i.*Wext*Dext(j)).*...
        FM_freq .* phiK(h*Wext));
end

FM_space = real(FM_space)./max(abs(real(FM_space)));

FM_space = align_to_reference(FM_space, f(Dext));
FM_space = FM_space(Dd-dd+1:Dd+dd)'; 
recspace = align_to_reference(recspace, f(D))';

l8_space_err = max(abs(recspace-f(D)))/max(abs(f(D)))
l8_alg_err = max(abs(FM_space-f(D)))/max(abs(f(D)))

%% FINAL PLOTTING 
figure
plot(Wext(length(Wext)/2:end), real(plotting(length(plotting)/2:end)))
hold on
plot(Wext(length(Wext)/2:end), real(recfft(length(recfft)/2:end)))
% plot(Wext(length(Wext)/2:end), imag(recfft(length(recfft)/2:end)))

figure
plot(D, recspace)
hold on
plot(D,FM_space)
hold on
plot(D,(f(D)))
legend('ours', 'FM', 'true')

toc
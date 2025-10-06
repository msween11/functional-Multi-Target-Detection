 rng(2);
symfft = @(x) ifftshift(fft(fftshift(x)));

N = 1; l = 5; B = 10;
x=-N:2^(-l):N-2^(-l);
X=-2*N:2^(-l):2*N-2^(-l);
w = -pi*(2^l):(pi/N):pi*(2^l)-(pi/N); 
W = -pi*(2^l):(pi/(2*N)):pi*(2^l)-(pi/(2*N)); 
Xext = -2*N*B:2^(-l):2*N*B-2^(-l);
Wext = -pi*(2^l):(pi/(2*N*B)):pi*(2^l)-(pi/(2*N*B)); 
WWext = -2*pi*(2^l):(pi/(2*N*B)):2*pi*(2^l)-(pi/(2*N*B)); 
wgap = Wext(2) - Wext(1);
L = length(Wext);

lambda = .1; 
sigma = .5; 
n = 100000;
h = 0.05;
r = 0.01;
zero_cond = 0;
thresh = 0.001;

% LINEAR SPLINE
% f1 = @(x) max(1-abs(x),0);
% f1 = @(x) f1(x) + f1(pi*x); 
% nrm = @(x) max(abs(f1(x)));
% f = @(x) f1(x) / nrm(x);
% f=@(x) f(x-.3);

% CUBIC SPLINE
     f2 = @(x) (1/6)*( 0.* (x<0)+...
    (x.^3) .* (0<= x & x < 1)+...
    (-3*(x-1).^3 + 3*(x-1).^2+3*(x-1)+1) .* (1<= x & x < 2)+...
    (3*(x-2).^3-6*(x-2).^2+4) .* (2 <= x & x < 3)+...
    ((4-x).^3) .* (3 <=x & x < 4)+...
    0 .* (x >= 4));
f2 = @(x) f2(2*x+2);
f2 = @(x) f2(x) + f2(pi*x);
nrm = @(x) max(abs(f2(x)));
f = @(x) f2(x) / nrm(x);

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



% UNIF SHIFTS
dist = makedist('Uniform', 'lower', -1, 'upper', 1);
shifts = random(dist,1, n);

% REJECTION SAMPLING
% mu = @(x) max(1-abs(x),0);
% mu = @(x) mu(x) + mu(pi*x);
% nrm_mu = @(x) sum(mu(x));
% mu = @(x) mu(x) / nrm_mu(x);
% 
% shifts = zeros(1,n);
% a = -1; b = 1; c= max(mu(x))+.1; p = 1; 
% while p < n+1
% 
%     U = rand; V = rand;
% 
%     Z = a+U*(b-a); Y = c*V;
% 
%     if Y<=mu(Z)
%        shifts(p)=Z;
%        p = p+1;
%     end   
% end

%% BEGIN RECOVERY NOW
% sampling in space

covr = zeros(length(X));
for j = 1:length(X)
    for k = 1:length(X)
        covr(j,k) = (sigma^2)*exp(-norm(X(j)-X(k))^2/(2*lambda^2)); 
    end
end

noise = mvnrnd(zeros(1,length(covr)), covr, n);

samples = zeros(n, length(X));
for j = 1:n
    samples(j,:) = f(X-shifts(j))+noise(j,:); 
end

n1 = length(X);
n2 = length(Xext);
nside = (n2-n1)/2;
covr_padded = [zeros(nside,n2); zeros(n1,nside) covr zeros(n1,nside); zeros(nside,n2)];
F = dftmtx(length(Xext));
cov_freq = center_BS(F'*center_BS(covr_padded)*F);

% constructing psihat in freq.
%%
samples_padded = [ zeros(n,(length(Xext)-length(samples(1,:)))/2)...
        samples...
        zeros(n,(length(Xext)-length(samples(1,:)))/2) ];

fsamples_padded = ifftshift(fft(fftshift(samples_padded,2),[],2),2);
psihat = conj((1/n*(fsamples_padded'*fsamples_padded))-cov_freq)/L;
clear('fsamples_padded')
% Doing Kotlarski's identity EMPIRICALLY

[empx, empy] = gradient(psihat, wgap);

adx = diag(empx); 
ad_emp = diag(psihat);

integrand_emp = adx./ad_emp; integrand_emp = integrand_emp';

rec_emp = zeros(1,length(Wext)/2);
rec_emp(1) = 0;
for j = 2:length(Wext)/2
    rec_emp(j) = rec_emp(j-1) + wgap*integrand_emp( length(Wext)/2 +j);
end

psireg = psihat./(L*min(ones(length(Wext))/L, r*sqrt(n)*abs(psihat)));

ad_reg = diag(psireg);
integrand_reg = adx./ad_reg; integrand_reg = integrand_reg';

rec_reg = zeros(1,length(Wext)/2);
rec_reg(1) = 0;
for j = 2:length(Wext)/2
    rec_reg(j) = rec_reg(j-1) + wgap*integrand_reg( length(Wext)/2 +j);
end

rec1_reg = (exp(rec_reg)); 

if(zero_cond == 1)
    growth_filt = .001;
    T = thresh;
    ps = real(ad_emp)';
    g = gradient(ps, wgap);
    gg = gradient(g, wgap);

    pidz = findz(g);
    pidz_filt = pidz(abs(gg(pidz))> growth_filt );
    pidz_filt = pidz_filt(Wext(pidz_filt)>0);
    pz_filt = Wext(pidz_filt);
    % fliploc = pz_filt(1:2:end); %this is by parity
    %^^^ all assumes that f^ft(0) = 1...else change even/odd
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

plotting = symfft(f(Xext));
plotting = plotting ./plotting(length(Xext)/2+1);
recfft = [conj(flip(rec1_reg)) (rec1_reg)];
l8_freq_rel_err = max(abs(recfft - plotting))/max(abs(plotting));
%hfix recovery in space now
phiK = @(w) 1 .* (-1 <= w & w <= 1);
%%

rec = zeros(1, length(Xext));
for j = 1:length(Xext)
    rec(j) = (wgap/(2*pi)).* sum(exp(1i.*Wext*Xext(j)).*...
        recfft .* phiK(h*Wext));
end
t = length(X)/2;  
tt = length(Xext)/2;
rec_trunc = real(rec(tt-t+1:tt+t)); %choosing real here should be ok
recspace = rec_trunc/max(abs(rec_trunc));

l8_hfix_space_err = max(abs(recspace-f(X)))/max(abs(f(X)));

% optimized h now
hrange = 1:.05:2;
h = 10.^(-hrange);

h_errs = zeros(1,length(h));

for w=h

idh = find(h == w);

rec = zeros(1, length(Xext));
for j = 1:length(Xext)
    rec(j) = (wgap/(2*pi)).* sum(exp(1i.*Wext*Xext(j)).*...
        recfft .* phiK(h(idh)*Wext));
end

rec_trunc_temp = real(rec(tt-t+1:tt+t)); %choosing real here should be ok
recspace_temp = rec_trunc_temp/max(abs(rec_trunc_temp));
h_errs(idh) = 100*max(abs(recspace_temp-f(X)))/max(abs(f(X)));

end

%---->h optimization over, now pick the optima
h_opt = h(find(h_errs == min(h_errs), 1, 'first'));
rec = zeros(1, length(Xext));
for j = 1:length(Xext)
    rec(j) = (wgap/(2*pi)).* sum(exp(1i.*Wext*Xext(j)).*...
        recfft .* phiK(h_opt*Wext));
end
rec_trunc = real(rec(tt-t+1:tt+t)); %choosing real here should be ok
recspace_opt = rec_trunc/max(abs(rec_trunc));

l8_hopt_space_err = max(abs(recspace_opt-f(X)))/max(abs(f(X)));

%% plotting 
figure
plot(Wext(length(Wext)/2:end), real(plotting(length(plotting)/2:end)))
hold on
plot(Wext(length(Wext)/2:end), real(recfft(length(recfft)/2:end)))
plot(Wext(length(Wext)/2:end), imag(recfft(length(recfft)/2:end)))

figure
plot(X, recspace_opt)
hold on
% plot(X, recspace)
plot(X, f(X))

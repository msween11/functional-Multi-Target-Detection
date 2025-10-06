rng(1);
symfft = @(x) ifftshift(fft(fftshift(x)));

% SIGNAL
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

%% ----REJECTION SAMPLING

n = 5; shifts = zeros(1,n);

for i = 2:n
    a = -3*i; b = 3*i; Z= [];
    while isempty(Z) == 1
        Z = a + (b-a).*rand;
        if min(abs(Z-shifts)) < 3
            Z = [];
        end
        Z(Z>=-3*(i-1) & Z<=3*(i-1)) = [] ; %enforcing density
    end
    shifts(i) = Z;
end
% shifts


%%  Constructing M(t)
sigma = .1; lambda = 0.1;
covr = @(x) (sigma^2)*exp(-abs(x-x').^2./(2*lambda^2)); 
noise = @(x) mvnrnd(zeros(1,length(covr(x))),covr(x), 1);

M = @(x) zeros(1,length(x));
for i = 1:n
    M = @(x) M(x) + f(x-shifts(i));
end

M = @(x) M(x)+ noise(x);

% MTD PLOT FOR PROPOSAL


mtdplot = tiledlayout(1,3, 'TileSpacing','tight')
ax1 = nexttile;

plot(ax1, D,f(D), LineWidth=2)
axis square
ax2 = nexttile;
plot(ax2, x, M(x), LineWidth=2)
axis square

sigma = 1; lambda = 0.1;
covr = @(x) (sigma^2)*exp(-abs(x-x').^2./(2*lambda^2)); 
noise = @(x) mvnrnd(zeros(1,length(covr(x))),covr(x), 1);

M = @(x) zeros(1,length(x));
for i = 1:n
    M = @(x) M(x) + f(x-shifts(i));
end

M = @(x) M(x)+ noise(x);
ax3 = nexttile;
plot(ax3, x,M(x), LineWidth=2)
axis square


%% Constructing A_3M(z_1,z_2)
N = 3*(n+1); l = 5; B = 10;
x=-N:2^(-l):N-2^(-l);
X=-2*N:2^(-l):2*N-2^(-l);
D = -3:2^(-l):3-2^(-l);
Dext = -3*B:2^(-l):3*B-2^(-l);
Wext = -pi*(2^l):(pi/(3*B)):pi*(2^l)-(pi/(3*B)); 
wgap = Wext(2) - Wext(1);
L = length(Wext);
%%
disp('data evaluation starting')
%evaluate M on extra large window X to make shifting easier.
%function evaluation is SLOW
datastrip = M(X);
disp('data evaluation done')
%shifting is just moving the window
data = zeros(length(x), length(D));
for i=1:length(D)
    j = D(i)*2^l;
    data(:,i) = datastrip(length(X)/4+j:3*length(X)/4-1+j);
end
%%
%lags / integrating range have to be carefully chosen....
A3M = zeros(length(D));

for i = 1:length(D)
    for  j = 1:length(D)
        A3M(i,j) = 2^(-l)*sum(data(:,length(D)/2+1).*data(:,i).*data(:,j));
    end
end
clear('data');

%%
A3M = (1/n)*A3M;

n1 = length(D);
n2 = length(Dext);
nside = (n2-n1)/2;
A3M_padded = [zeros(nside,n2); zeros(n1,nside) A3M zeros(n1,nside); zeros(nside,n2)];

F = dftmtx(length(Dext));
A3M_ft = center_BS(F'*center_BS(A3M_padded)*F);

%% Unbiasing for A3M
u = integral(f,-1,1)*sigma^2*n;
u = 1/u;

%% DOING KOTLARSKI NOW

psihat = A3M_ft;
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

r = 1; h = 0.02;
psireg = psihat./(L*min(ones(length(Wext))/L, r*sqrt(n)*abs(psihat)));

ad_reg = diag(psireg);
integrand_reg = adx./ad_reg; integrand_reg = integrand_reg';

rec_reg = zeros(1,length(Wext)/2);
rec_reg(1) = 0;
for j = 2:length(Wext)/2
    rec_reg(j) = rec_reg(j-1) + wgap*integrand_reg( length(Wext)/2 +j);
end

rec1_reg = (exp(rec_reg)); 

plotting = symfft(f(Dext));
plotting = plotting ./plotting(length(Dext)/2+1);
recfft = [conj(flip(rec1_reg)) (rec1_reg)];
l8_freq_rel_err = max(abs(recfft - plotting))/max(abs(plotting));

%hfix recovery in space now
phiK = @(w) 1 .* (-1 <= w & w <= 1);


rec = zeros(1, length(Dext));
for j = 1:length(Dext)
    rec(j) = (wgap/(2*pi)).* sum(exp(1i.*Wext*Dext(j)).*...
        recfft .* phiK(h*Wext));
end
t = length(D)/2;  
tt = length(Dext)/2;
rec_trunc = real(rec(tt-t+1:tt+t)); %choosing real here should be ok
recspace = rec_trunc/max(abs(rec_trunc));

l8_hfix_space_err = max(abs(recspace-f(D)))/max(abs(f(D)));

%% FINAL PLOTTING 
figure
plot(Wext(length(Wext)/2:end), real(plotting(length(plotting)/2:end)))
hold on
plot(Wext(length(Wext)/2:end), real(recfft(length(recfft)/2:end)))
plot(Wext(length(Wext)/2:end), imag(recfft(length(recfft)/2:end)))

figure
hold on
plot(D, recspace)
plot(D, f(D))


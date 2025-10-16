rng(1);
n = 1000; 
N = 3*(n+1); l = 5; B = 10;
x=-N:2^(-l):N-2^(-l);
X=-2*N:2^(-l):2*N-2^(-l);
D = -3:2^(-l):3-2^(-l);
Dext = -3*B:2^(-l):3*B-2^(-l);
Wext = -pi*(2^l):(pi/(3*B)):pi*(2^l)-(pi/(3*B)); 
wgap = Wext(2) - Wext(1);
L = length(Wext);

symfft = @(x) ifftshift(fft(fftshift(x)));

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
%% SAMPLING SHIFTS

s = -3 + 6*rand(1, n); s(1) = abs(s(1)); %enforcing parity
shifts = zeros(1,n);

evenidx = 2:2:n;
shifts(evenidx) = -3*evenidx + s(evenidx - 1);

oddidx = 3:2:n;
shifts(oddidx) = 3*oddidx + s(oddidx - 1);

%%  CONSTRUCTING M(t)...slowest chunk...not sure how to speed up?
M = zeros(1,length(X));

for j = 1:n 
    Y = X - shifts(j);
    [~, c] = min(abs(Y));
    M = M + [zeros(1,c-64) f(Y(c-64:c+64)) ...
        zeros(1,length(X)-length(zeros(1,c-64)) - length(f(Y(c-64:c+64))))];
    % disp(j)
end


%% CONSTRUCTING NOISE 

sigma = .5; lambda = 1;
rho = @(h) (sigma^2)*exp(-abs(h(1)-h(2)).^2./(2*lambda^2)); 
ep = stationary_Gaussian_process(1,length(X), rho); 
datastrip = M + ep;

%% CONSTRUCTING A3M
%shifting is just moving the window
data = zeros(length(x), length(D));
for i=1:length(D)
    j = D(i)*2^l;
    data(:,i) = datastrip(length(X)/4+j:3*length(X)/4-1+j);
end

v = data(:, length(D)/2 + 1);   
A3M = (1/n)*2^(-l) * ((data .* v)' * data); 

%% UNBIASING... DOES NOT WORK
% int = integral(f,-1,1);
% bias = sigma^2*int; 
% covr = @(x) (sigma^2)*exp(-abs(x).^2./(2*lambda^2)); 
% temp = zeros(1,length(D));
% for i = 1:length(D)
%     temp(i) = rho([D(length(D)/2+1),D(i)]);
% end
% temp2 = covr(D).*integral(f,-1,1);
% 
% for i = 1:length(D)
%     A3M(i,:) = A3M(i,:)-(integral(f,-1,1)*temp);
% end
% 
% A3M(1:(length(D)+1):length(D)^2) =...
%     A3M(1:(length(D)+1):length(D)^2) - bias;
% 
% A3M(:,length(D)/2+1) = A3M(:,length(D)/2+1) - bias;
% A3M(length(D)/2+1,:) = A3M(length(D)/2+1,:) - bias;
% %check larger lambda, shouldn't be constant......
% figure
% imagesc(A3M)
% 


%% PADDING A3M AND MAPPING INTO FREQUENCY 

n1 = length(D);
n2 = length(Dext);
nside = (n2-n1)/2;
A3M_padded = [zeros(nside,n2); zeros(n1,nside) A3M zeros(n1,nside); zeros(nside,n2)];

F = dftmtx(length(Dext));
A3M_ft = center_BS(F'*center_BS(A3M_padded)*F);


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

r = .0001; h = 0.05;
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


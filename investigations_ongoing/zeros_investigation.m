
rng(1);

[params, grids] = setparams(n= 2^10, l = 5, B = 10, z = 1, ...
    t = .01, sigma = .0001, lambda = .1, r = .1, h = .05);

covr = @(x) (params.sigma^2)*exp(-abs(x).^2./(2*params.lambda^2)); 
covr_freq = real(grids.dx*symfft(covr(grids.Dext)));
%both covr and covr_freq are scaled correctly according to cont theory

func_choice = 'f4';
f = choose_function(func_choice); 

if strcmp(func_choice, 'f4')
   params.zero_cond = 1;
else 
   params.zero_cond = 0;
end

idx = find(grids.Dext == 0);
func = f(grids.Dext);

fft = real(symfft(func));
phases = angle(fft(idx:end));
phase_deriv = gradient(phases);

% figure
% plot(grids.Dext(idx:end), fft(idx:end))
% hold on
% plot(grids.Dext(idx:end),5*phases)
% legend('fft', 'phase')


%so it seems like if we can get a good approximation for
% phase_deriv, then we should be able to recover the phases and thus flips,
% so let's see if we can recover phase_deriv!

%%

%getting our observation strip
datastrip = get_obs(params, grids, f, type = 'discrete');

%constructing our invariants
[~, psihat,psireg] = get_psihat_auto(params, grids, f, datastrip);

%% manual frequency marching
B = center_BS(psireg);
N = size(B, 1);
B = (B+B')/2;
z0 = sign(B(1, 1));
z1 = 1;

Psi = angle(B);

psi = zeros(floor((N+1)/2), 1);
psi(1) = angle(z0);
psi(2) = angle(z1);

for j = 3 : ceil((N+1)/2)

    psi_est = zeros(ceil(j/2)-1,1); %blank vector to be averaged over

    for k=2:ceil(j/2) %multiple paths
        psi_est(k-1) =  Psi( j, k ) + psi(k) + psi(j-k+1); %FM recursion
    end

    % Averaging all estimation over SO(2).
    psi(j) = angle(sum(exp(1j*psi_est)));

    %OR, taking max
    % [~, idx] = max(abs(psi_est)); 
    % psi(j) = psi_est(idx);
end

if N/2 == floor(N/2)
    psi = [angle(z0) ; psi(2:end)  ; flipud(-psi(2:end-1))];
else
    psi = [angle(z0) ; psi(2:end)  ; flipud(-psi(2:end))];    
end

zout = exp(1i*psi);

%now our code
dd = grids.lenD/2;
Dd = length(grids.Dext)/2;
p = circshift(real(diag(psireg)), Dd); %choose psihat here maybe?

FM_freq = sqrt(p).*zout;
FM_freq = circshift(FM_freq', Dd);



%% kotlarski below here

%to make this path variable, we will need to make this a loop over freq and
%replace the cumsum with a sum
mid = find(grids.Wext == 0);
Wpos = grids.Wext(mid:end);

%initialization
kot = zeros(1,length(Wpos));

%need to make paths and thus npaths

for j = 1:length(Wpos)
    kot_est = zeros(1,npaths);
    for k = 2:ceil(j/2)
        % psi_est(k-1) =  Psi( j, k ) + psi(k) + psi(j-k+1); %FM recursion
        
    end

    kot(j) = angle(sum(exp(1j*kot_est)));
end


empx = gradient(psihat, grids.dw);

adx = diag(empx); %this is the path choice, change this
ad_emp = diag(psihat); %and here

integrand_emp = adx./ad_emp; integrand_emp = integrand_emp';

f_emp = integrand_emp(mid:end);
rec_emp = cumsum(f_emp)*grids.dw;

ad_reg = diag(psireg);
integrand_reg = adx./ad_reg; integrand_reg = integrand_reg';

f_reg = integrand_reg(mid:end);
rec_reg = cumsum(f_reg)*grids.dw;
rec1_reg = (exp(rec_reg)); 
%%

test_pd = imag(f_reg);
test_p = cumsum(test_pd)*grids.dw;

% figure; hold on;
% plot(Wpos, test_pd);
% plot(Wpos, fft(idx:end))



% figure
% plot(test_p)
% 



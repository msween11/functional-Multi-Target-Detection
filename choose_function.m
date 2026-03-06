function f = choose_function(signal)
    
if(signal == 'f1')
f1 = @(x) max(1-abs(x),0);
f1 = @(x) f1(x) + f1(pi*x); 
nrm = @(x) max(abs(f1(x)));
f = @(x) f1(x) / nrm(x);

elseif(signal == 'f2')
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

elseif(signal == 'f3')
k=30;
f3 = @(x) exp(-20*x.^2).*(cos(k*x)+.5);
% f3 = @(x) f3(x-.3);
nrm = @(x) max(abs(f3(x)));
f = @(x) f3(x) / nrm(x);

elseif(signal == 'f4')
means = [-.4;.4]; devs = [.02];
f4 = gmdistribution(means, devs);
f4 = @(x) pdf(f4, x')';
nrm = @(x) max(abs(f4(x)));
f = @(x) f4(x) / nrm(x);
% f = @(x) f(x-.1);

elseif(signal == 'f5')
k=8;
f3 = @(x) exp(-20*x.^2).*(cos(k*x));
f3 = @(x) f3(x-.3);
nrm = @(x) max(abs(f3(x)));
f = @(x) f3(x) / nrm(x);

end



function [A1M, A2M,A3M, psihat, psireg] = get_psihat(params, grids, f, datastrip)


d = length(grids.X)-length(grids.x); d=d/2;
datastrip_trunc = datastrip(d:end-d-1);

data = zeros(length(grids.x), length(grids.D)); %data for A3M
    
for k = 1:grids.lenD
        j = grids.D(k)/grids.dx;
        data(:,k) = datastrip(d+j:end-d-1+j);
end
   
A1M =  (1/params.n)*grids.dx * sum(datastrip_trunc);
A2M =  (1/params.n)*grids.dx * (datastrip_trunc * data);   
A3M = (1/params.n)*grids.dx * ((data .* datastrip_trunc')' * data); 

%  UNBIASING  

covr = @(x) (params.sigma^2)*exp(-abs(x).^2./(2*params.lambda^2)); 
u = covr(grids.D)*integral(f,-1,1);
A3M = A3M-u-u'-(integral(f,-1,1)*covr(grids.D-grids.D'));

% PADDING A3M AND MAPPING INTO FREQUENCY 

n1 = grids.lenD;
n2 = length(grids.Dext);
nside = (n2-n1)/2;
A3M_padded = [zeros(nside,n2); zeros(n1,nside) A3M zeros(n1,nside); zeros(nside,n2)];

F = dftmtx(n2);
A3M_ft = center_BS(F'*center_BS(A3M_padded)*F);
psihat = A3M_ft./A1M.*grids.dx^2;
nrm = 1/(A1M.*grids.dx^2);

psireg = psihat./(nrm*min(ones(length(grids.Wext))/nrm,...
    params.r*sqrt(params.n)*abs(psihat)));


end
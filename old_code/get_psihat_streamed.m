function [A3M, psihat, psireg] = get_psihat_streamed(params, grids, f, datastrip)

d = length(grids.X)-length(grids.x); d=d/2;
datastrip_trunc = datastrip(d:end-d-1);

% Compute lower triangle over full D x D lag grid, then reflect once.
% Avoids non-contiguous upper-triangle writes inside the inner loop.
A3M = zeros(grids.lenD, grids.lenD);

for i = 1:grids.lenD
    ii = grids.D(i)/grids.dx;
    di = datastrip_trunc .* datastrip(d+ii:end-d-1+ii);
    for j = 1:i
        jj = grids.D(j)/grids.dx;
        A3M(i,j) = 1/params.n*grids.dx*dot(di, datastrip(d+jj:end-d-1+jj));
    end
end
A3M = A3M + tril(A3M,-1)';

%  UNBIASING
covr = @(x) (params.sigma^2)*exp(-abs(x).^2./(2*params.lambda^2));
u    = covr(grids.D)*integral(f,-1,1);
A3M  = A3M - u - u' - (integral(f,-1,1)*covr(grids.D-grids.D'));

% PADDING A3M AND MAPPING INTO FREQUENCY
n1    = grids.lenD;
n2    = length(grids.Dext);
nside = (n2-n1)/2;
A3M_padded = [zeros(nside,n2); zeros(n1,nside) A3M zeros(n1,nside); zeros(nside,n2)];

F      = dftmtx(n2);
A3M_ft = center_BS(F'*center_BS(A3M_padded)*F);

A1M    = (1/params.n)*grids.dx * sum(datastrip_trunc);
psihat = A3M_ft ./ A1M .* grids.dx^2;

nrm    = 1/(A1M.*grids.dx^2);
psireg = psihat ./ (nrm*min(ones(length(grids.Wext))/nrm, ...
    params.r*sqrt(params.n)*abs(psihat)));

end

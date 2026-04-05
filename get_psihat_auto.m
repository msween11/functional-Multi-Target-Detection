function [ps, psihat, psireg, psihat_gradx] = get_psihat_auto(params, grids, f, datastrip)
% Combines get_psihat and get_psihat_streamed, choosing between them based
% on available system memory checked BEFORE allocation.
%
% On HPC Linux systems the OS kills the process before MATLAB can throw an
% out-of-memory error, so try-catch is unreliable. Instead we query
% /proc/meminfo (Linux) or the MATLAB memory() function (Windows) up front
% and pick the streamed path if the data matrix would not fit safely.

d = length(grids.X)-length(grids.x); d=d/2;
datastrip_trunc = datastrip(d:end-d-1);

A1M = (1/params.n)*grids.dx * sum(datastrip_trunc);

% Memory needed for the full data matrix (lenx x lenD, float64)
data_bytes = length(datastrip_trunc) * grids.lenD * 8;

if enough_memory(data_bytes)

    % ---- NON-STREAMED: build data matrix once, vectorised A2M and A3M ----
    data = zeros(length(grids.x), grids.lenD);
    for k = 1:grids.lenD
        j = grids.D(k)/grids.dx;
        data(:,k) = datastrip(d+j:end-d-1+j);
    end

    A2M = (1/params.n)*grids.dx * (datastrip_trunc * data);
    A3M = (1/params.n)*grids.dx * ((data .* datastrip_trunc')' * data);

else

    % ---- STREAMED: lower-triangle dot-product loop, no large matrix ----
    fprintf('get_psihat_auto: insufficient memory for data matrix, using streamed path.\n');

    % A2M: cheap loop, no full matrix needed
    A2M = zeros(1, grids.lenD);
    for k = 1:grids.lenD
        j = grids.D(k)/grids.dx;
        A2M(k) = (1/params.n)*grids.dx * dot(datastrip_trunc, datastrip(d+j:end-d-1+j));
    end

    % A3M: compute lower triangle, reflect
    A3M = zeros(grids.lenD, grids.lenD);
    fprintf('A3M (streamed):   0%%');
    for i = 1:grids.lenD
        fprintf('\b\b\b\b%3d%%', ...
            round(100 * i*(i+1) / (grids.lenD*(grids.lenD+1))));
        ii = grids.D(i)/grids.dx;
        di = datastrip_trunc .* datastrip(d+ii:end-d-1+ii);
        for j = 1:i
            jj = grids.D(j)/grids.dx;
            A3M(i,j) = 1/params.n*grids.dx*dot(di, datastrip(d+jj:end-d-1+jj));
        end
    end
    fprintf('\n');
    A3M = A3M + tril(A3M,-1)';

end

% ---- SHARED: unbiasing, padding, frequency map, psihat ----

covr = @(x) (params.sigma^2)*exp(-abs(x).^2./(2*params.lambda^2));
covr_freq = real(grids.dx*symfft(covr(grids.Dext)));

u    = covr(grids.D)*integral(f,-1,1);
A3M  = A3M - u - u' - (integral(f,-1,1)*covr(grids.D-grids.D'));

n1    = grids.lenD;
n2    = length(grids.Dext);
nside = (n2-n1)/2;
A3M_padded = [zeros(nside,n2); zeros(n1,nside) A3M zeros(n1,nside); zeros(nside,n2)];
A2M_padded = [zeros(1,nside) A2M zeros(1,nside)];

F      = dftmtx(n2);
A3M_ft = center_BS(F'*center_BS(A3M_padded)*F);
psihat = A3M_ft./A1M.*grids.dx^2;

% Horizontal gradient of psihat via FT derivative property.
% The row direction in center_BS(F'*center_BS(*)*F) uses a positive
% exponent e^{+i*w1*z1}, so differentiation w.r.t. w1 brings down +i*z1:
%   d/dw1 psihat(w1,w2) = FT[+i*z1 * A3(z1,z2)](w1,w2) / A1M * dx^2
% Multiply A3M row-wise by +i*z1 (rows index z1 = grids.D),
% then apply the same 2D transform and normalisation.
B3M          = (1i * grids.D(:)) .* A3M;
B3M_padded   = [zeros(nside,n2); zeros(n1,nside) B3M zeros(n1,nside); zeros(nside,n2)];
psihat_gradx = center_BS(F'*center_BS(B3M_padded)*F) ./ A1M .* grids.dx^2;
nrm    = 1/(A1M.*grids.dx^2);

psireg = psihat./(nrm*min(ones(length(grids.Wext))/nrm,...
    params.r*sqrt(params.n)*abs(psihat)));

constant = (max(grids.X)-min(grids.X))/params.n;
ps = real(grids.dx*symfft(A2M_padded))-constant*covr_freq;


end


function ok = enough_memory(required_bytes)
% Returns true if required_bytes is safely within available system memory.
% Reads /proc/meminfo on Linux (HPC), falls back to MATLAB memory() on Windows.
    safety = 0.75;  % only use up to 75% of what's available
    try
        if isunix
            [~, txt] = system('grep MemAvailable /proc/meminfo');
            available_kb = sscanf(txt, 'MemAvailable: %d');
            available_bytes = available_kb * 1024;
        else
            [~, sys] = memory();
            available_bytes = sys.PhysicalMemory.Available;
        end
        ok = required_bytes < safety * available_bytes;
    catch
        % If memory query fails, fall back to streamed to be safe
        ok = false;
    end
end

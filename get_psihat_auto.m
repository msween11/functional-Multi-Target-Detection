function [A1M, A2M, A3M, psihat, psireg] = get_psihat_auto(params, grids, f, datastrip)
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
    for i = 1:grids.lenD
        ii = grids.D(i)/grids.dx;
        di = datastrip_trunc .* datastrip(d+ii:end-d-1+ii);
        for j = 1:i
            jj = grids.D(j)/grids.dx;
            A3M(i,j) = 1/params.n*grids.dx*dot(di, datastrip(d+jj:end-d-1+jj));
        end
    end
    A3M = A3M + tril(A3M,-1)';

end

% ---- SHARED: unbiasing, padding, frequency map, psihat ----

covr = @(x) (params.sigma^2)*exp(-abs(x).^2./(2*params.lambda^2));
u    = covr(grids.D)*integral(f,-1,1);
A3M  = A3M - u - u' - (integral(f,-1,1)*covr(grids.D-grids.D'));

n1    = grids.lenD;
n2    = length(grids.Dext);
nside = (n2-n1)/2;
A3M_padded = [zeros(nside,n2); zeros(n1,nside) A3M zeros(n1,nside); zeros(nside,n2)];

F      = dftmtx(n2);
A3M_ft = center_BS(F'*center_BS(A3M_padded)*F);
psihat = A3M_ft./A1M.*grids.dx^2;
nrm    = 1/(A1M.*grids.dx^2);

psireg = psihat./(nrm*min(ones(length(grids.Wext))/nrm,...
    params.r*sqrt(params.n)*abs(psihat)));

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

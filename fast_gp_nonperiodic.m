function gp = fast_gp_nonperiodic(Nx, dx, sigma, lambda)
    % Nx: number of grid points
    % dx: spatial step size
    % sigma, lambda: covariance parameters
    % Output gp has length Nx

    % simulate on a 2x domain to avoid periodic wrap-around
    Nbig = 2 * Nx;
    freq = ifftshift((-Nbig/2:Nbig/2-1) / (Nbig*dx));
    omega = 2*pi*freq;

    % Squared exponential spectral density
    S = sigma^2 * sqrt(2*pi) * lambda .* exp(-0.5 * (lambda^2) .* (omega.^2));

    % Hermitian-symmetric random spectrum
    half = Nbig/2;
    F = zeros(1,Nbig);
    pos = 1:(half+1);
    neg = (half+2):Nbig;

    wp = (randn(1,length(pos)) + 1i*randn(1,length(pos))) / sqrt(2);
    F(pos) = sqrt(S(pos)) .* wp;
    F(neg) = conj(F(pos(end-1:-1:2)));

    % Inverse FFT to real space
    gp_big = real(ifft(F)) * sqrt(Nbig);

    % Crop the center to avoid edge artifacts
    start_idx = Nbig/4 + 1;
    end_idx   = start_idx + Nx - 1;
    gp = gp_big(start_idx:end_idx);

    % Normalize variance
    gp = gp * (sigma / std(gp));
end
function [mean_est, P_est, B_est] = invariants_from_data(X, sigma)
% Given M noisy, circularly shifted observations of a signal of length N in
% a matrix X of size N x M and the standard deviation sigma of the noise,
% returns an estimate of the mean, the power spectrum and the bispectrum of
% the signal. Uses a parfor for the bispectrum.
%
% For the power spectrum, the noise power is subtracted to avoid bias. If
% this results in negative entries, those are trimmed to zero.
%
% For the bispectrum, it is an estimate of the bispectrum of the centered
% signal (mean removed) to avoid bias. This operation affects only the
% first row and column of B as well as its diagonal.
%
% If sigma is not specified, it will be estimated from the data.
%
% May 2017
% https://arxiv.org/abs/1705.00641
% https://github.com/NicolasBoumal/MRA

    [N, M] = size(X);
    
    %% If sigma is not supplied, estimate it from the data
    if ~exist('sigma', 'var') || isempty(sigma)
        sigma = std(sum(X, 1))/sqrt(N);
    end

    %% Estimate the mean and subtract it from the observations.
    %  This also gives an estimate of the first component of the fft of x.
    %  Centering the observations helps for the bispectrum estimation part.
    mean_est = mean(X(:));
    Xc = X - mean_est;

    %% Prepare DFTs of centered signals
    Xc_fft = fft(Xc);

    %% Estimate the power spectrum (gives estimate of modulus of DFT of x).
    P_est = mean(abs(Xc_fft).^2, 2);
    % Debias the estimate using the supplied or estimated value of sigma
    P_est = P_est - N*sigma^2;
    % Take the nonnegative part
    P_est = max(0, P_est);
    
    %% Estimate the bispectrum
    if nargout >= 3
        
        B_est = zeros(N);
        parfor m = 1 : M
            xm_fft = Xc_fft(:, m);
            Bm = (xm_fft*xm_fft') .* circulant(xm_fft);
            B_est = B_est + Bm;
        end
        B_est = B_est/M;
        
    end

end

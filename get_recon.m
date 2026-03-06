function X_recon = get_recon(Emean,power,phase)
        X_recon = sqrt(power).*phase; 
        X_recon = real(ifft(X_recon));
%         X_recon = align_chen(X_recon,Emean);
        X_recon = X_recon - mean(X_recon)  + Emean;
end
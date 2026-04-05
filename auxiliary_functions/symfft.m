function out = symfft(x)
    out = ifftshift(fft(fftshift(x)));
end

function [X, Y, E] = Gaussian_Beam_ASM(z_distance, N, L)
    lambda = 1;         
    k = 2*pi/lambda;    
    
    dx = L / N;
    x = (-N/2 : N/2 - 1) * dx;
    y = x;
    [X, Y] = meshgrid(x, y);

    w0 = 1.6 * lambda;
    %w0 = sqrt(2) / (0.1 * k);
    E_z0 = exp(-(X.^2 + Y.^2) / (w0^2));
    

    dk = 2*pi/L;
    kx_v = (-N/2 : N/2 - 1) * dk;
    ky_v = kx_v;
    [Kx, Ky] = meshgrid(kx_v, ky_v);

    Kz_sq = k^2 - Kx.^2 - Ky.^2;
    Kz = sqrt(Kz_sq);
    Kz(Kz_sq < 0) = -1i * sqrt(abs(Kz_sq(Kz_sq < 0)));

    Propagator = exp(-1i * Kz * z_distance);

    A_z0 = fftshift(fft2(E_z0));
    
    A_z = A_z0 .* Propagator;
    
    E = ifft2(ifftshift(A_z));

end
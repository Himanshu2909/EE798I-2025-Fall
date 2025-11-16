function [X, Y, E] = Gaussian_Beam_Matrix(z, N, L)
lambda = 1;
k = 2*pi/lambda;
sigma = 0.1*k;

dk = 2*pi/L;
kx = (-N/2 : N/2 - 1) * dk;
ky = kx;
[Kx, Ky] = meshgrid(kx, ky);

%{
Issue #1
Kz = (k^2 - Kx.^2 - Ky.^2).^0.5;

Does not handle the case k^2 < Kx.^2 - Ky.^2 well (evanescent waves)
For this case, sqrt() or .^0.5 returns positive imaginary root, hence
exp(-1i * Kz * z) becomes exp(-1i * (Ki) * z) = exp(Kz) which is 
physically incorrect as the wave can't be exponentially growing. 

%}
Kz_square = k^2 - Kx.^2 - Ky.^2;
Kz = sqrt(Kz_square);
Kz(Kz_square < 0) = -1i * sqrt(abs(Kz_square(Kz_square < 0)));

A = exp(-(Kx.^2 + Ky.^2)/(2*sigma^2));

dx = L / N;
x = (-N/2 : N/2 - 1) * dx;
y = x;
[X, Y] = meshgrid(x, y);

%{
for i1 = 1 : size(Kz,1)
    for i2 =  1 : size(Kz,1)
        Atmp = A(i1,i2)*exp(-1i*Kz(i1,i2)*z*lambda/450);
        E = E + Atmp*exp(-1i*(Kx(i1,i2).*X + Ky(i1,i2).*Y));
    end
end

Nested for loop method is computationally inefficient, with a complexity of
N(K_z^2)
%}

Kx_v = Kx(:);
Ky_v = Ky(:); 

Propagator = exp(-1i * Kz * z);

% sum is scaled by dkx * dky, giving better approximation
dkx = kx(2) - kx(1);
dky = ky(2) - ky(1);
A_z = A .* Propagator * (dkx * dky);
A_z_v = A_z(:); 

X_v = X(:); 
Y_v = Y(:);

Phase_Matrix = exp(-1i * (X_v * Kx_v.' + Y_v * Ky_v.'));
E_v = Phase_Matrix * A_z_v;

E = reshape(E_v, size(X));

w0 = sqrt(2) / (0.1 * k);
normalization_factor = (w0^2) / (4 * pi);
E = E * normalization_factor;

end


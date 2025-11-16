function [X, Y, E_tmp] = Gaussian_Beam(z)

lambda = 1; k = 2*pi/lambda; sigma = 0.1*k;

kx = linspace(-0.3, 0.3, 50)*k;
ky = linspace(-0.3, 0.3, 50)*k;
[Kx, Ky] = meshgrid(kx, ky);
Kz = (k^2 - Kx.^2 - Ky.^2).^0.5;
A = exp(-(Kx.^2 + Ky.^2)/(2*sigma^2));

x = linspace(-3, 3, 30)*lambda;
y = linspace(-3, 3, 30)*lambda;
[X, Y] = meshgrid(x, y);

E = zeros(size(X));

for i1 = 1 : size(Kz,1)
    for i2 =  1 : size(Kz,1)
        Atmp = A(i1,i2)*exp(-1i*Kz(i1,i2)*z*lambda/450);
        E = E + Atmp*exp(-1i*(Kx(i1,i2).*X + Ky(i1,i2).*Y));
    end
end

E_tmp = real(E);

end


z_distance = 20;
lambda = 1;
N = 64; % Number of grid points ( must be power of 2 for FFT )
L = 30 * lambda ;

[X, Y, E_complex] = Gaussian_Beam(z_distance);
%[X, Y, E_complex] = Gaussian_Beam_Matrix(z_distance, N, L);

Intensity = abs(E_complex).^2;
Phase = angle(E_complex);

figure;

% Plot 1: 2D Intensity (Heatmap)
subplot(1, 3, 1);
pcolor(X, Y, Intensity);
shading interp;           
axis equal tight;         
title(['Intensity at z = ', num2str(z_distance)]);
xlabel('x / \lambda');
ylabel('y / \lambda');
colorbar;

% Plot 2: 3D Intensity (Surface)
subplot(1, 3, 2);
surf(X, Y, Intensity);
shading interp;           
title(['Intensity at z = ', num2str(z_distance)]);
xlabel('x / \lambda');
ylabel('y / \lambda');
zlabel('Intensity');

% Plot 3: 2D Phase (Wavefront)
subplot(1, 3, 3);
pcolor(X, Y, Phase);
shading interp;
axis equal tight;
title(['Phase at z = ', num2str(z_distance)]);
xlabel('x / \lambda');
ylabel('y / \lambda');
colorbar;
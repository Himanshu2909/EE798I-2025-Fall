function plot_poynting_vector()
    z_distance = 20;    
    lambda = 1;        
    k = 2*pi/lambda;    
    
    N = 256;           
    L = 30 * lambda;
    dx = L / N;
    dy = dx;
    
    dz = lambda / 20;
    
    [X, Y, E_x] = Gaussian_Beam_ASM(z_distance, N, L);
    
    [~, ~, E_x_plus] = Gaussian_Beam_ASM(z_distance + dz, N, L);
    
    [~, ~, E_x_minus] = Gaussian_Beam_ASM(z_distance - dz, N, L);
    
    x_v = (-N/2 : N/2 - 1) * dx;
    y_v = x_v;

    [~, dEx_dy] = gradient(E_x, dx, dy);
    
    dEx_dz = (E_x_plus - E_x_minus) / (2 * dz);
    
    H_y = (1 / (-1i * k)) * dEx_dz;
    H_z = (1 / (-1i * k)) * (-dEx_dy);
    
    
    S_x = zeros(size(E_x)); 
    S_y = 0.5 * real(-E_x .* conj(H_z));
    S_z = 0.5 * real(E_x .* conj(H_y));
    
    figure('Name', 'Poynting Vector Component Magnitudes', 'Position', [100, 100, 1000, 400]);
    
    subplot(1, 2, 1);
    imagesc(x_v, y_v, S_z);
    axis image;
    colorbar;
    title(sprintf('S_z (Propagation) at z = %.1f', z_distance), 'FontSize', 12);
    xlabel('x / \lambda');
    ylabel('y / \lambda');
    
    subplot(1, 2, 2);
    imagesc(x_v, y_v, S_y);
    axis image;
    colorbar;
    title(sprintf('S_y (Transverse) at z = %.1f', z_distance), 'FontSize', 12);
    xlabel('x / \lambda');
    ylabel('y / \lambda');
    
    max_S_z = max(S_z, [], 'all');
    max_S_y = max(abs(S_y), [], 'all');
    fprintf('Max S_z (propagation): %.4e\n', max_S_z);
    fprintf('Max S_y (transverse):  %.4e\n', max_S_y);
    fprintf('Ratio (S_z / S_y):     %.2f\n', max_S_z / max_S_y);
    
    
    step = 16;
    idx = 1:step:N;
    X_q = X(idx, idx);
    Y_q = Y(idx, idx);
    S_x_q = S_x(idx, idx); 
    S_y_q = S_y(idx, idx);
    
    S_mag = sqrt(S_x_q.^2 + S_y_q.^2);
    S_mag(S_mag == 0) = 1; 
    
    figure('Name', 'Transverse Poynting Vector (S_x, S_y)', 'Position', [200, 200, 700, 600]);
    imagesc(x_v, y_v, S_z);
    hold on;
    
    quiver(X_q, Y_q, S_x_q./S_mag, S_y_q./S_mag, 'r', 'LineWidth', 1);
    
    hold off;
    axis image;
    colorbar;
    title('Transverse Power Flow (S_x, S_y) overlaid on S_z', 'FontSize', 12);
    xlabel('x / \lambda');
    ylabel('y / \lambda');
    legend('Transverse S-Vector Direction');
end
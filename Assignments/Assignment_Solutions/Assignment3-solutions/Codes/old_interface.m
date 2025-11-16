function old_interface(P)
    n1 = P.n1;
    n2 = P.n2;
    z_interface = P.z_interface;
    z_obs = P.z_obs;
    lambda_vac = P.lambda_vac;
    w0 = P.w0;
    theta_inc_deg = P.theta_inc_deg;
    polarization = P.polarization;
    N = P.N;
    L = P.L;
    case_title = P.case_title;

    k0 = 2*pi / lambda_vac;  
    k1 = n1 * k0;            
    k2 = n2 * k0;            
    theta_inc_rad = deg2rad(theta_inc_deg);
    
    dx = L / N;
    x = (-N/2 : N/2 - 1) * dx;
    y = x;
    [X, Y] = meshgrid(x, y);
    
    dk = 2*pi / L;
    kx_v = (-N/2 : N/2 - 1) * dk;
    ky_v = kx_v;
    [Kx, Ky] = meshgrid(kx_v, ky_v);
    
    E_z0 = exp(-(X.^2 + Y.^2) / (w0^2));
    E_z0_tilted = E_z0 .* exp(-1i * k1 * sin(theta_inc_rad) * X);
    A_z0_tilted = fftshift(fft2(E_z0_tilted));
    
    Kz1_sq = k1^2 - Kx.^2 - Ky.^2;
    Kz2_sq = k2^2 - Kx.^2 - Ky.^2;
    
    Kz1 = sqrt(Kz1_sq);
    Kz1(Kz1_sq < 0) = -1i * sqrt(abs(Kz1_sq(Kz1_sq < 0)));
    
    Kz2 = sqrt(Kz2_sq);
    Kz2(Kz2_sq < 0) = -1i * sqrt(abs(Kz2_sq(Kz2_sq < 0)));
    
    Cos_t1 = Kz1 / k1;
    Cos_t2 = Kz2 / k2;
    
    epsilon = 1e-9;
    
    % *** BUG FIX ***
    % Added definition for den_s, which was missing in the original code.
    den_s = (n1 .* Cos_t1 + n2 .* Cos_t2) + epsilon; 
    den_p = (n2 .* Cos_t1 + n1 .* Cos_t2) + epsilon;
    
    if strcmpi(polarization, 's')
        R = (n1 .* Cos_t1 - n2 .* Cos_t2) ./ den_s;
        T = (2 .* n1 .* Cos_t1) ./ den_s;
    elseif strcmpi(polarization, 'p')
        R = (n2 .* Cos_t1 - n1 .* Cos_t2) ./ den_p;
        T = (2 .* n1 .* Cos_t1) ./ den_p;
    else
        error("Polarization must be 's' or 'p'.");
    end
    
    if z_obs < z_interface
        % Calculate total field (Incident + Reflected)
        P_inc = exp(-1i * Kz1 * z_obs);
        A_inc_at_z_obs = A_z0_tilted .* P_inc;
        
        P_to_iface = exp(-1i * Kz1 * z_interface);
        A_refl_at_iface = A_z0_tilted .* P_to_iface .* R;
        
        P_back_from_iface = exp(+1i * Kz1 * (z_obs - z_interface)); 
        A_refl_at_z_obs = A_refl_at_iface .* P_back_from_iface;
        
        A_final = A_inc_at_z_obs + A_refl_at_z_obs;
        
        plot_title_prefix = 'Total Field (Incident + Reflected)';
        
    else
        % Calculate transmitted field
        P_to_iface = exp(-1i * Kz1 * z_interface);
        A_trans_at_iface = A_z0_tilted .* P_to_iface .* T;
        
        dz_2 = z_obs - z_interface;
        P_in_med2 = exp(-1i * Kz2 * dz_2);
        
        A_final = A_trans_at_iface .* P_in_med2;
        
        plot_title_prefix = 'Transmitted Field';
    end
    
    E_final = ifft2(ifftshift(A_final));
    
    Intensity = abs(E_final).^2;
    Phase = angle(E_final);
    
    % Get y=0 cross-section (ensure correct indexing for even N)
    Intensity_cross_section = Intensity(N/2 + 1, :); 
    
    % --- Plotting ---
    figure('Name', case_title, 'Position', [100, 100, 1300, 500]);
    
    % Add an overall title to the figure window
    sgtitle(case_title, 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    subplot(1, 3, 1);
    imagesc(x, y, Intensity);
    axis image;
    colorbar;
    title(sprintf('%s Intensity at z = %.1f', plot_title_prefix, z_obs), 'FontSize', 12);
    xlabel('x / \lambda_{vac}');
    ylabel('y / \lambda_{vac}');
    if z_obs < z_interface
        text(x(5), y(end-10), sprintf('...interface at z=%.1f', z_interface), 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    subplot(1, 3, 2);
    imagesc(x, y, Phase);
    axis image;
    colorbar;
    title(sprintf('%s Phase at z = %.1f', plot_title_prefix, z_obs), 'FontSize', 12);
    xlabel('x / \lambda_{vac}');
    ylabel('y / \lambda_{vac}');
    
    subplot(1, 3, 3);
    plot(x, Intensity_cross_section, 'LineWidth', 2);
    title(sprintf('Intensity Cross-Section (y=0) at z = %.1f', z_obs), 'FontSize', 12);
    xlabel('x / \lambda_{vac}');
    ylabel('Intensity (a.u.)');
    grid on;
    xlim([min(x) max(x)]);
    
    fprintf('Simulation complete for: %s\n', case_title);
end
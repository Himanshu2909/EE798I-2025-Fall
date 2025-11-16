function interface_cases()
    % Main function to define and run all simulation cases.
    
    % --- Define Default Parameters ---
    % These will be used as a base for all cases
    P_base.n1 = 1.0;            
    P_base.n2 = 1.5;           
    P_base.z_interface = 15;    
    P_base.z_obs = 14;          
    P_base.lambda_vac = 1.0;    
    P_base.w0 = 3 * P_base.lambda_vac; 
    P_base.theta_inc_deg = 20;  
    P_base.polarization = 'p';  
    P_base.N = 256;             
    P_base.L = 50 * P_base.lambda_vac; 
    P_base.case_title = 'Default Case';

    fprintf('Running all simulation cases...\n\n');
    
    % --- Case 1: Standard Reflection & Interference ---
    P_case1 = P_base;
    P_case1.case_title = 'Case 1: Standard Reflection (z < z_{iface})';
    P_case1.n1 = 1.0;
    P_case1.n2 = 1.5;
    P_case1.theta_inc_deg = 20;
    P_case1.z_obs = 10; % Observe reflection
    P_case1.polarization = 's';
    plot_beam_at_interface(P_case1);

    % --- Case 2: Standard Transmission ---
    P_case2 = P_base;
    P_case2.case_title = 'Case 2: Standard Transmission (z > z_{iface})';
    P_case2.n1 = 1.0;
    P_case2.n2 = 1.5;
    P_case2.theta_inc_deg = 20;
    P_case2.z_obs = 20; % Observe transmission
    plot_beam_at_interface(P_case2);
    
    % --- Case 3: Brewster's Angle (p-pol, no reflection) ---
    P_case3 = P_base;
    theta_b = atand(P_case3.n2 / P_case3.n1); % Calculate Brewster's angle
    P_case3.case_title = sprintf('Case 3: Brewster Angle (%.1f°) - p-pol (No Reflection)', theta_b);
    P_case3.theta_inc_deg = theta_b;
    P_case3.polarization = 'p';
    P_case3.z_obs = 10; % Observe reflection (should be none)
    plot_beam_at_interface(P_case3);

    % --- Case 4: Brewster's Angle (s-pol, for comparison) ---
    P_case4 = P_case3; % Start from case 3 parameters
    P_case4.case_title = sprintf('Case 4: Brewster Angle (%.1f°) - s-pol (Strong Reflection)', theta_b);
    P_case4.polarization = 's'; % Just change polarization
    plot_beam_at_interface(P_case4);

    % --- Case 5: Total Internal Reflection (TIR) ---
    P_case5 = P_base;
    P_case5.n1 = 1.5; % From dense
    P_case5.n2 = 1.0; % To sparse
    theta_c = asind(P_case5.n2 / P_case5.n1); % Critical angle
    P_case5.theta_inc_deg = theta_c + 5; % Go 5 degrees above critical angle
    P_case5.case_title = sprintf('Case 5: TIR (\\theta_c=%.1f°) - Reflected Beam', theta_c);
    P_case5.z_obs = 10; % Observe reflected beam
    plot_beam_at_interface(P_case5);

    % --- Case 6: Evanescent Wave (TIR) ---
    P_case6 = P_case5; % Start from TIR case
    P_case6.case_title = sprintf('Case 6: TIR (\\theta_c=%.1f°) - Evanescent Wave', theta_c);
    P_case6.z_obs = 15.2; % Observe JUST past interface
    plot_beam_at_interface(P_case6);

    % --- Case 7: Evanescent Wave Decay (TIR) ---
    P_case7 = P_case5; % Start from TIR case
    P_case7.case_title = sprintf('Case 7: TIR (\\theta_c=%.1f°) - Evanescent Decay', theta_c);
    P_case7.z_obs = 17.0; % Observe further past interface
    plot_beam_at_interface(P_case7);
    
    fprintf('\nAll simulations complete.\n');
end
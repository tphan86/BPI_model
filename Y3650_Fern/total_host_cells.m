function T = total_host_cells(para_fit,t)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pre-estimate parameters for the model
    K = 10^9; % CFU/mL
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fitting parameters for the model
    r = para_fit(1);
    a = para_fit(2);
    a_ss = para_fit(3);
    a_rr = para_fit(4);
    delta = para_fit(5);
    gamma = para_fit(6);
    phi = para_fit(7);
    beta = para_fit(8);
    m = para_fit(9);
    S0 = para_fit(10);
    R0 = para_fit(11);
    P0 = para_fit(12);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the total host cells from the model
    initial_value = [S0;R0;P0];
    p = [r,a,a_ss,a_rr,K,delta,gamma,phi,beta,m];
    func = @one_species_one_phage;
    % Set ODE solver options to handle stiffness and improve stability
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [~,y] = ode15s(@(r,y) func(r,y,p),t,initial_value,options);
    T = y(:,1) + y(:,2);
end
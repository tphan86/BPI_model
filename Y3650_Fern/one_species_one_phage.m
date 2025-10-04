function dy = one_species_one_phage(t,y,p)
    %%%%%%%%%%%%%%%
    r = p(1); % growth rate of S in the absence of P and R
    a = p(2); % interspecific competition coefficient of S and R
    a_ss = p(3); % intraspecific competition coefficient of S
    a_rr = p(4); % intraspecific competition coefficient of R
    K = p(5); % system-wide carrying capacity
    delta = p(6); % rate of resistance acquisition
    gamma = p(7); % fitness cost of phage resistance (rel. to parental phenotype)
    phi = p(8); % the adsorption rate of P when attaching to S
    beta = p(9); % burst size of phage P when infecting S
    m = p(10); % decay rate of phage P
    %%%%%%%%%%%%%%%
    S = y(1); % the density of phage-susceptible bacteria S
    R = y(2); % the density of phage-resistant bacteria R
    P = y(3); % the density of phage P
    %%%%%%%%%%%%%%%
    dy = zeros(3,1);
    dy(1) = r*S*(1-(a_ss*S + a*R)/K) - delta*S - phi*P*S;
    dy(2) = gamma*r*R*(1-(a*S + a_rr*R)/K) + delta*S;
    dy(3) = beta*phi*P*S - m*P;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
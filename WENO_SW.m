function [rho, u, p, e] = WENO_SW(dx, rho0, u0, p0, tEnd)
%input initial condition & total time
global gamma CFL N epsilon
t = 0;
rho = rho0;
u = u0;
p = p0;
E = p ./ ((gamma-1) * rho) + u.^2/2;
U = [rho; rho.*u; rho.*E];   % U values at j
while t <= tEnd
    % update time step
    a = sqrt(gamma * p ./ rho);
    S_max = max(max(abs(u) + a));
    dt = 0.0001;
    % time advance
    U_sub = U;
    % 1st step
    L = Calculate_L_WENO_SW(U, dx);
    U = U_sub + dt * L;
    % 2nd step
    L = Calculate_L_WENO_SW(U, dx);
    U = 3/4 * U_sub + 1/4 * (U + dt * L);
    % 3rd step
    L = Calculate_L_WENO_SW(U, dx);
    U =1/3 * U_sub + 2/3 * (U + dt * L);
    % compute flow properties
    rho = U(1,:);
    u = U(2,:)./rho;
    E = U(3,:)./rho;
    p = (gamma-1).*(E-u.^2./2).*rho;
    e = p./(gamma-1)./rho;
    t = t+dt;
end
end
function  L = Calculate_L_WENO_SW(U, dx)
    global gamma CFL N epsilon
    rho = U(1,:);
    u = U(2,:)./rho;
    E = U(3,:)./rho;
    p = (gamma-1).*(E-u.^2./2).*rho;
    H = (p * gamma) ./ ((gamma-1) * rho) + u.^2/2;
    a = sqrt(gamma * p ./ rho);
    % compute eigenvalues
    lambda_1 = u - a;
    lambda_2 = u;
    lambda_3 = u + a;
    % split the flux
    lambda_pos_1 = (lambda_1 + abs(lambda_1)) / 2;
    lambda_pos_2 = (lambda_2 + abs(lambda_2)) / 2;
    lambda_pos_3 = (lambda_3 + abs(lambda_3)) / 2;
    lambda_neg_1 = (lambda_1 - abs(lambda_1)) / 2;
    lambda_neg_2 = (lambda_2 - abs(lambda_2)) / 2;
    lambda_neg_3 = (lambda_3 - abs(lambda_3)) / 2;

    lambda_pos_1 = repmat(lambda_pos_1,3,1);
    lambda_pos_2 = repmat(lambda_pos_2,3,1);
    lambda_pos_3 = repmat(lambda_pos_3,3,1);
    lambda_neg_1 = repmat(lambda_neg_1,3,1);
    lambda_neg_2 = repmat(lambda_neg_2,3,1);
    lambda_neg_3 = repmat(lambda_neg_3,3,1);
    
    K_1 = [ones(1,N+3); u-a; H-u.*a];
    K_2 = [ones(1,N+3); u; 1/2 * u.^2];
    K_3 = [ones(1,N+3); u+a; H+u.*a];
    
    F_pos = lambda_pos_1 .* rho / 2 / gamma .* K_1 + lambda_pos_2 .* rho * (gamma-1) / gamma .* K_2 + lambda_pos_3 .* rho / 2 / gamma .* K_3;
    F_neg = lambda_neg_1 .* rho / 2 / gamma .* K_1 + lambda_neg_2 .* rho * (gamma-1) / gamma .* K_2 + lambda_neg_3 .* rho / 2 / gamma .* K_3;
    
    % WENO
    % j+1/2
    IS_pos_1 = 1/4 * (F_pos(:,1:N+1) - 4*F_pos(:,2:N+2)+ 3*F_pos(:,3:N+3)).^2 + 13/12 * (F_pos(:,1:N+1) - 2*F_pos(:,2:N+2) + F_pos(:,3:N+3)).^2;
    IS_pos_2 = 1/4 * (F_pos(:,1:N+1) - F_pos(:,3:N+3)).^2 + 13/12 * (F_pos(:,1:N+1) - 2*F_pos(:,2:N+2) + F_pos(:,3:N+3)).^2;
    IS_pos_3 = 1/4 * (3*F_pos(:,1:N+1) - 4*F_pos(:,2:N+2)+ F_pos(:,3:N+3)).^2 + 13/12 * (F_pos(:,1:N+1) - 2*F_pos(:,2:N+2) + F_pos(:,3:N+3)).^2;
    
    alpha_pos_1 = 1/10./ (ones(3,N+1) * 10^(-6) + IS_pos_1).^2; 
    alpha_pos_1 = [zeros(3,1) zeros(3,1) alpha_pos_1];
    alpha_pos_2 = 6/10./ (ones(3,N+1) * 10^(-6) + IS_pos_2).^2; 
    alpha_pos_2 = [zeros(3,1) alpha_pos_2 zeros(3,1)];
    alpha_pos_3 = 3/10./ (ones(3,N+1) * 10^(-6) + IS_pos_3).^2; 
    alpha_pos_3 = [alpha_pos_3 zeros(3,1) zeros(3,1)];
    
    omega_pos_1 = alpha_pos_1 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);
    omega_pos_2 = alpha_pos_2 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);
    omega_pos_3 = alpha_pos_3 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);
    
    F_pos_1 = 1/3 * F_pos(:,1:N+1) - 7/6 * F_pos(:,2:N+2) + 11/6 * F_pos(:,3:N+3);
    F_pos_1 = [zeros(3,1) zeros(3,1) F_pos_1];
    F_pos_2 = -1/6 * F_pos(:,1:N+1) + 5/6 * F_pos(:,2:N+2) + 1/3 * F_pos(:,3:N+3);
    F_pos_2 = [zeros(3,1) F_pos_2 zeros(3,1)];
    F_pos_3 = 1/3 * F_pos(:,1:N+1) + 5/6 * F_pos(:,2:N+2) - 1/6 * F_pos(:,3:N+3);
    F_pos_3 = [F_pos_3 zeros(3,1) zeros(3,1)];
    
    F_pos_c = omega_pos_1 .* F_pos_1 + omega_pos_2 .* F_pos_2 + omega_pos_3 .* F_pos_3;
    
    % j-1/2
    IS_neg_1 = 1/4 * (F_neg(:,3:N+3) - 4*F_neg(:,2:N+2) + 3*F_neg(:,1:N+1)).^2 + 13/12 * (F_neg(:,3:N+3) - 2*F_neg(:,2:N+2) + F_neg(:,1:N+1)).^2;
    IS_neg_2 = 1/4 * (F_neg(:,3:N+3) - F_neg(:,1:N+1)).^2 + 13/12 * (F_neg(:,3:N+3) - 2*F_neg(:,2:N+2) + F_neg(:,1:N+1)).^2;
    IS_neg_3 = 1/4 * (3*F_neg(:,3:N+3) - 4*F_neg(:,2:N+2)+ F_neg(:,1:N+1)).^2 + 13/12 * (F_neg(:,3:N+3) - 2*F_neg(:,2:N+2) + F_neg(:,1:N+1)).^2;
    
    alpha_neg_1 = 1/10./ (ones(3,N+1) * 10^(-6) + IS_neg_1).^2; 
    alpha_neg_1 = [alpha_neg_1 zeros(3,1) zeros(3,1)];
    alpha_neg_2 = 6/10./ (ones(3,N+1) * 10^(-6) + IS_neg_2).^2; 
    alpha_neg_2 = [zeros(3,1) alpha_neg_2 zeros(3,1)];
    alpha_neg_3 = 3/10./ (ones(3,N+1) * 10^(-6) + IS_neg_3).^2; 
    alpha_neg_3 = [zeros(3,1) zeros(3,1) alpha_neg_3]; 

    
    omega_neg_1 = alpha_neg_1 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);
    omega_neg_2 = alpha_neg_2 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);
    omega_neg_3 = alpha_neg_3 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);

    F_neg_1 = 1/3 * F_neg(:,3:N+3) - 7/6 * F_neg(:,2:N+2) + 11/6 * F_neg(:,1:N+1);
    F_neg_1 = [F_neg_1 zeros(3,1) zeros(3,1)];
    F_neg_2 = -1/6 * F_neg(:,3:N+3) + 5/6 * F_neg(:,2:N+2) + 1/3 * F_neg(:,1:N+1);
    F_neg_2 = [zeros(3,1) F_neg_2 zeros(3,1)];
    F_neg_3 = 1/3 * F_neg(:,3:N+3) + 5/6 * F_neg(:,2:N+2) - 1/6 * F_neg(:,1:N+1);
    F_neg_3 = [zeros(3,1) zeros(3,1) F_neg_3];
    
    F_neg_c = omega_neg_1 .* F_neg_1 + omega_neg_2 .* F_neg_2 + omega_neg_3 .* F_neg_3;
    
    for i = 2:N+2
        L(:,i) = - 1/dx * (F_pos_c(:,i) - F_pos_c(:,i-1) + F_neg_c(:,i+1) - F_neg_c(:,i));
    end
    L(:,1) = 0;      %transmissive left boundary
    L(:,N+3) = 0;  %transmissive right boundary
end
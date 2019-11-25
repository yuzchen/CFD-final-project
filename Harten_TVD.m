function [rho, u, p, e] = Harten_TVD(dx, rho0, u0, p0, tEnd, entropy_fix)
%input initial condition & total time
global gamma CFL N epsilon
t = 0;
rho = rho0;
u = u0;
p = p0;
E = p ./ ((gamma-1) * rho) + u.^2/2;
U = [rho; rho.*u; rho.*E];   % U values at i
limiter_type = input('please input limiter type: \n {1} minmod; {2} superbee.\n');
% limiter_type = 1;
while t <= tEnd
    % update time step using CFL condition
    a = sqrt(gamma * p ./ rho);
    S_max = max(max(abs(u) + a));
    dt = dx * CFL / S_max;
    % time advance using RK3
    % 1st step
    U_sub = U;
    L = Calculate_L_TVD(U, dt, dx, limiter_type);
    U = U_sub + dt * L;
    % 2nd step
    L = Calculate_L_TVD(U, dt, dx, limiter_type);
    U = 3/4 * U_sub + 1/4 * (U + dt * L);
    % 3rd step
    L = Calculate_L_TVD(U, dt, dx, limiter_type);
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
function L = Calculate_L_TVD(U, dt, dx, limiter_type)
    global gamma CFL N epsilon
    rho = U(1,:);
    u = U(2,:)./rho;
    E = U(3,:)./rho;
    p = (gamma-1).*(E-u.^2./2).*rho;
    H = (p * gamma) ./ ((gamma-1) * rho) + u.^2/2;
    a = sqrt(gamma * p ./ rho);
    % Roe average
    rho_bar = (sqrt(rho(1:N+2)) + sqrt(rho(2:N+3)))/2;
    u_bar = (sqrt(rho(1:N+2)) .* u(1:N+2) + sqrt(rho(2:N+3)) .* u(2:N+3)) ./ (sqrt(rho(1:N+2)) + sqrt(rho(2:N+3)));
    H_bar = (sqrt(rho(1:N+2)) .* H(1:N+2) + sqrt(rho(2:N+3)) .* H(2:N+3)) ./ (sqrt(rho(1:N+2)) + sqrt(rho(2:N+3)));
    a_bar = sqrt((gamma-1) * (H_bar - u_bar.^2 / 2));
    % averaged eigenvalues
    lambda_bar_1 = u_bar - a_bar;
    lambda_bar_2 = u_bar;
    lambda_bar_3 = u_bar + a_bar;
    % averaged left and right eigenvectors
    r_bar_1 = [ones(1,N+2); u_bar-a_bar; H_bar-u_bar.*a_bar];
    r_bar_2 = [ones(1,N+2); u_bar.^2; 1/2 * u_bar.^2];
    r_bar_3 = [ones(1,N+2); u_bar+a_bar; H_bar+u_bar.*a_bar];
    l_bar_1 = [u_bar ./ a_bar /4 .* (2*ones(1,N+2) + (gamma-1) * u_bar./a_bar); -[ones(1,N+2) + (gamma-1)*u_bar./a_bar]./a_bar/2; ones(1,N+2)*(gamma-1)/2./a_bar.^2];
    l_bar_2 = [ones(1,N+2) - (gamma-1)/2 * u_bar.^2 ./ a_bar.^2; (gamma-1)*u_bar.^2 ./ a_bar.^2; -ones(1,N+2)*(gamma-1)./a_bar.^2];
    l_bar_3 = [-u_bar ./ a_bar /4 .* (2*ones(1,N+2) - (gamma-1) * u_bar./a_bar); [ones(1,N+2) - (gamma-1)*u_bar./a_bar]./a_bar/2; ones(1,N+2)*(gamma-1)/2./a_bar.^2];
    % calculate alpha_{j+1/2} and \tilde{g}_{j+1/2}
    alpha_bar_1 = l_bar_1 .* (U(:,2:N+3) - U(:,1:N+2));
    alpha_bar_1 = sum(alpha_bar_1,1);
    alpha_bar_2 = l_bar_2 .* (U(:,2:N+3) - U(:,1:N+2));
    alpha_bar_2 = sum(alpha_bar_2,1);
    alpha_bar_3 = l_bar_3 .* (U(:,2:N+3) - U(:,1:N+2));
    alpha_bar_3 = sum(alpha_bar_3,1);
    g_bar_1 = 1/2 * (Q(dt/dx*lambda_bar_1, epsilon) - (dt/dx*lambda_bar_1).^2) .* alpha_bar_1; 
    g_bar_2 = 1/2 * (Q(dt/dx*lambda_bar_2, epsilon) - (dt/dx*lambda_bar_2).^2) .* alpha_bar_2; 
    g_bar_3 = 1/2 * (Q(dt/dx*lambda_bar_3, epsilon) - (dt/dx*lambda_bar_3).^2) .* alpha_bar_3; 
    % apply limiter
    g_1 = zeros(1,N+3);
    g_2 = zeros(1,N+3);
    g_3 = zeros(1,N+3);
    switch limiter_type
        case{1} % minmod
            for i = 2:N+2
                g_1(i) = minmod([g_bar_1(i) g_bar_1(i-1)]);
                g_2(i) = minmod([g_bar_2(i) g_bar_2(i-1)]);
                g_3(i) = minmod([g_bar_3(i) g_bar_3(i-1)]);
            end
        case{2} % superbee
            for i = 2:N+2
                g_1(i) = minmod([2*g_bar_1(i) 2*g_bar_1(i-1) (g_bar_1(i)+g_bar_1(i-1))/2]);
                g_2(i) = minmod([2*g_bar_2(i) 2*g_bar_2(i-1) (g_bar_2(i)+g_bar_2(i-1))/2]);
                g_3(i) = minmod([2*g_bar_3(i) 2*g_bar_3(i-1) (g_bar_3(i)+g_bar_3(i-1))/2]);
            end   
    otherwise
        error('wrong limiter type');
    end
    for i = 1:N+2
        if alpha_bar_1(i) >= 1e-6
            gamma_bar_1(i) = (g_1(i+1)-g_1(i)) / alpha_bar_1(i);
            gamma_bar_2(i) = (g_2(i+1)-g_2(i)) / alpha_bar_2(i);
            gamma_bar_3(i) = (g_3(i+1)-g_3(i)) / alpha_bar_3(i);
        else
            gamma_bar_1(i) = 0;
            gamma_bar_2(i) = 0;
            gamma_bar_3(i) = 0;
        end
    end
    % Psi_{j+1/2}
    Psi_bar_1 = dx/dt * (g_1(2:N+3)-g_1(1:N+2) - Q(dt/dx*lambda_bar_1 + gamma_bar_1, epsilon) .* alpha_bar_1);
    Psi_bar_2 = dx/dt * (g_2(2:N+3)-g_2(1:N+2) - Q(dt/dx*lambda_bar_2 + gamma_bar_2, epsilon) .* alpha_bar_2);
    Psi_bar_3 = dx/dt * (g_3(2:N+3)-g_3(1:N+2) - Q(dt/dx*lambda_bar_3 + gamma_bar_3, epsilon) .* alpha_bar_3);
    Psi_bar_1 =  repmat(Psi_bar_1,3,1);
    Psi_bar_2 =  repmat(Psi_bar_2,3,1);
    Psi_bar_3 =  repmat(Psi_bar_3,3,1);
    % numerical flux
    F = [rho.*u; rho.*u.^2+p; rho.*u.*H]; % F(U_j^n)
    F_c = 1/2 * (F(:,1:N+2) + F(:,2:N+3) + Psi_bar_1 .* r_bar_1 + Psi_bar_2 .* r_bar_2 + Psi_bar_3 .* r_bar_3);
    % calculate L
    for i = 2:N+2
        L(:,i) = 1/dx * (F_c(:,i-1) - F_c(:,i));
    end
    L(:,1) = 0;      %transmissive left boundary
    L(:,N+3) = 0;  %transmissive right boundary
end
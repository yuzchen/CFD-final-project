function [rho, u, p, e] = SCM(dx, rho0, u0, p0, tEnd)
%input initial condition & total time
global gamma CFL N epsilon
t = 0;
rho = rho0;
u = u0;
p = p0;
E = p ./ ((gamma-1) * rho) + u.^2/2;
U = [rho; rho.*u; rho.*E];   % U values at j
while t <= tEnd
    % update time step using CFL condition
    a = sqrt(gamma * p ./ rho);
    S_max = max(max(abs(u) + a));
    dt = dx * CFL / S_max;
    % time advance
    U_sub = U;
    L = Calculate_L_SCM(U, dx);
    U = U_sub + dt * L;
%     L = Calculate_L_SCM(U, dx);
%     U = 3/4 * U_sub + 1/4 * (U + dt * L);
%     L = Calculate_L_SCM(U, dx);
%     U = 1/3 * U_sub + 2/3 * (U + dt * L);
    
    
    % compute flow properties
    rho = U(1,:);
    u = U(2,:)./rho;
    E = U(3,:)./rho;
    p = (gamma-1).*(E-u.^2./2).*rho;
    e = p./(gamma-1)./rho;
    t = t+dt;
end
end
function  L = Calculate_L_SCM(U, dx)
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
    lambda_bar_1 = u_bar - a_bar;
    lambda_bar_2 = u_bar;
    lambda_bar_3 = u_bar + a_bar;
    % lambda_{j+1/2}
    lambda_pos_1 = (lambda_bar_1 + abs(lambda_bar_1)) / 2;
    lambda_pos_2 = (lambda_bar_2 + abs(lambda_bar_2)) / 2;
    lambda_pos_3 = (lambda_bar_3 + abs(lambda_bar_3)) / 2;
    lambda_neg_1 = (lambda_bar_1 - abs(lambda_bar_1)) / 2;
    lambda_neg_2 = (lambda_bar_2 - abs(lambda_bar_2)) / 2;
    lambda_neg_3 = (lambda_bar_3 - abs(lambda_bar_3)) / 2;
    for i = 2:N+2
        % S_{j+1/2}^{-1}, S_{j+1/2}
        S_inv = [1                           1                      1;
                u_bar(i)-a_bar(i)            u_bar(i)               u_bar(i)+a_bar(i);
                H_bar(i)-u_bar(i)*a_bar(i)   1/2*u_bar(i)^2         H_bar(i)+u_bar(i)*a_bar(i)];
        S = [u_bar(i)/4/a_bar(i)*(2+(gamma-1)*u_bar(i)/a_bar(i))        -1/2/a_bar(i)*(1+(gamma-1)*u_bar(i)/a_bar(i))       (gamma-1)/2/a_bar(i)^2;
            1-(gamma-1)/2*u_bar(i)^2/a_bar(i)^2                         (gamma-1)*u_bar(i)/a_bar(i)^2                       -(gamma-1)/a_bar(i)^2;
            -u_bar(i)/4/a_bar(i)*(2-(gamma-1)*u_bar(i)/a_bar(i))        1/2/a_bar(i)*(1-(gamma-1)*u_bar(i)/a_bar(i))        (gamma-1)/2/a_bar(i)^2];
        % or just
%         S = inv(S_inv);
        lambda_pos = [lambda_pos_1(i)   0                     0;
                      0                 lambda_pos_2(i)       0;
                      0                 0                     lambda_pos_3(i)];
        lambda_neg = [lambda_neg_1(i)   0                     0;
                      0                 lambda_neg_2(i)       0;
                      0                 0                     lambda_neg_3(i)];
        for k = i-1:i+1
            V_pos(:,k) = lambda_pos * S * U(:,k);
            V_neg(:,k) = lambda_neg * S * U(:,k);
        end
        % calculate flux using upwind scheme
        V_pos_c = V_pos(:,i);   % V_{j+1/2}^+
        V_neg_c = V_neg(:,i+1); % V_{j+1/2}^-
        F(:,i) = S_inv * (V_neg_c + V_pos_c);
    end
    F(:,1) = F(:,2);
    for i = 2:N+2
        L(:,i) = 1/dx * (F(:,i-1) - F(:,i));
    end
    L(:,1) = 0;      %transmissive left boundary
    L(:,N+3) = 0;  %transmissive right boundary
end
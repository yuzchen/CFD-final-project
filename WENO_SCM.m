function [rho, u, p, e] = WENO_SCM(dx, rho0, u0, p0, tEnd)
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
%     S_max = max(max(abs(u) + a));
    dt = 0.0001;
    % time advance
    U_sub = U;
    % 1st step
    L = Calculate_L_WENO_SCM(U, dx);
    U = U_sub + dt * L;
    % 2nd step
    L = Calculate_L_WENO_SCM(U, dx);
    U = 3/4 * U_sub + 1/4 * (U + dt * L);
    % 3rd step
    L = Calculate_L_WENO_SCM(U, dx);
    U = 1/3 * U_sub + 2/3 * (U + dt * L);
    % compute flow properties
    rho = U(1,:);
    u = U(2,:)./rho;
    E = U(3,:)./rho;
    p = (gamma-1).*(E-u.^2./2).*rho;
    e = p./(gamma-1)./rho;
    t = t+dt;
end
end
function  L = Calculate_L_WENO_SCM(U, dx)
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
    V_pos = zeros(3,N+3);  V_neg = zeros(3,N+3);
    V_pos_c = zeros(3,N+2); V_neg_c = zeros(3,N+2);
    for i = 1:N+2
        % S_{j+1/2}^{-1}, S_{j+1/2}
%         S_inv = [1                           1                      1;
%                 u_bar(i)-a_bar(i)            u_bar(i)               u_bar(i)+a_bar(i);
%                 H_bar(i)-u_bar(i)*a_bar(i)   1/2*u_bar(i)^2         H_bar(i)+u_bar(i)*a_bar(i)];
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
        if i ~= 1 && i ~= 2 && i ~= N+2
            for k = i-2:i+2
                V_pos(:,k) = lambda_pos * S * U(:,k) ;   
                V_neg(:,k) = lambda_neg * S * U(:,k) ;
            end
            % WENO
            % j+1/2
            IS_pos_1 = 1/4 * (V_pos(:,i-2) - 4*V_pos(:,i-1)+ 3*V_pos(:,i)).^2 + 13/12 * (V_pos(:,i-2) - 2*V_pos(:,i-1) + V_pos(:,i)).^2;
            IS_pos_2 = 1/4 * (V_pos(:,i-1) - V_pos(:,i+1)).^2 + 13/12 * (V_pos(:,i-1) - 2*V_pos(:,i) + V_pos(:,i+1)).^2;
            IS_pos_3 = 1/4 * (3*V_pos(:,i) - 4*V_pos(:,i+1)+ V_pos(:,i+2)).^2 + 13/12 * (V_pos(:,i) - 2*V_pos(:,i+1) + V_pos(:,i+2)).^2;

            alpha_pos_1 = 1/10./ (ones(3,1) * epsilon + IS_pos_1).^2; 
            alpha_pos_2 = 6/10./ (ones(3,1) * epsilon + IS_pos_2).^2; 
            alpha_pos_3 = 3/10./ (ones(3,1) * epsilon + IS_pos_3).^2; 

            omega_pos_1 = alpha_pos_1 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);
            omega_pos_2 = alpha_pos_2 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);
            omega_pos_3 = alpha_pos_3 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);

            V_pos_1 = 1/3 * V_pos(:,i-2) - 7/6 * V_pos(:,i-1) + 11/6 * V_pos(:,i);
            V_pos_2 = -1/6 * V_pos(:,i-1) + 5/6 * V_pos(:,i) + 1/3 * V_pos(:,i+1);
            V_pos_3 = 1/3 * V_pos(:,i) + 5/6 * V_pos(:,i+1) - 1/6 * V_pos(:,i+2);
            % j-1/2
            IS_neg_1 = 1/4 * (V_neg(:,i+2) - 4*V_neg(:,i+1) + 3*V_neg(:,i)).^2 + 13/12 * (V_neg(:,i+2) - 2*V_neg(:,i+1) + V_neg(:,i)).^2;
            IS_neg_2 = 1/4 * (V_neg(:,i+1) - V_neg(:,i-1)).^2 + 13/12 * (V_neg(:,i+1) - 2*V_neg(:,i) + V_neg(:,i-1)).^2;
            IS_neg_3 = 1/4 * (3*V_neg(:,i) - 4*V_neg(:,i-1)+ V_neg(:,i-2)).^2 + 13/12 * (V_neg(:,i) - 2*V_neg(:,i-1) + V_neg(:,i-2)).^2;

            alpha_neg_1 = 1/10./ (ones(3,1) * epsilon + IS_neg_1).^2; 
            alpha_neg_2 = 6/10./ (ones(3,1) * epsilon + IS_neg_2).^2; 
            alpha_neg_3 = 3/10./ (ones(3,1) * epsilon + IS_neg_3).^2; 

            omega_neg_1 = alpha_neg_1 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);
            omega_neg_2 = alpha_neg_2 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);
            omega_neg_3 = alpha_neg_3 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);

            V_neg_1 = 1/3 * V_neg(:,i+2) - 7/6 * V_neg(:,i+1) + 11/6 * V_neg(:,i);
            V_neg_2 = -1/6 * V_neg(:,i+1) + 5/6 * V_neg(:,i) + 1/3 * V_neg(:,i-1);
            V_neg_3 = 1/3 * V_neg(:,i) + 5/6 * V_neg(:,i-1) - 1/6 * V_neg(:,i-2);
        elseif i == 1
            for k = i:i+2
                V_pos(:,k) = lambda_pos * S * U(:,k) ;
                V_neg(:,k) = lambda_neg * S * U(:,k) ;
            end
            % WENO
            % j+1/2
            IS_pos_3(:,i) = 1/4 * (3*V_pos(:,i) - 4*V_pos(:,i+1)+ V_pos(:,i+2)).^2 + 13/12 * (V_pos(:,i) - 2*V_pos(:,i+1) + V_pos(:,i+2)).^2;

            alpha_pos_1 = 0; 
            alpha_pos_2 = 0; 
            alpha_pos_3 = 3/10./ (ones(3,1) * epsilon + IS_pos_3).^2; 

            omega_pos_1 = alpha_pos_1 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);
            omega_pos_2 = alpha_pos_2 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);
            omega_pos_3 = alpha_pos_3 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);

            V_pos_1 = 0;
            V_pos_2 = 0;
            V_pos_3 = 1/3 * V_pos(:,i) + 5/6 * V_pos(:,i+1) - 1/6 * V_pos(:,i+2);
            % j-1/2
            IS_neg_1 = 1/4 * (V_neg(:,i+2) - 4*V_neg(:,i+1) + 3*V_neg(:,i)).^2 + 13/12 * (V_neg(:,i+2) - 2*V_neg(:,i+1) + V_neg(:,i)).^2;

            alpha_neg_1 = 1/10./ (ones(3,1) * epsilon + IS_neg_1).^2; 
            alpha_neg_2 = 0; 
            alpha_neg_3 = 0; 

            omega_neg_1 = alpha_neg_1 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);
            omega_neg_2 = alpha_neg_2 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);
            omega_neg_3 = alpha_neg_3 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);

            V_neg_1 = 1/3 * V_neg(:,i+2) - 7/6 * V_neg(:,i+1) + 11/6 * V_neg(:,i);
            V_neg_2 = 0;
            V_neg_3 = 0;
        elseif i == 2
            for k = i-1:i+2
                V_pos(:,k) = lambda_pos * S * U(:,k) ;
                V_neg(:,k) = lambda_neg * S * U(:,k) ;
            end
            % WENO
            % j+1/2
            IS_pos_2 = 1/4 * (V_pos(:,i-1) - V_pos(:,i+1)).^2 + 13/12 * (V_pos(:,i-1) - 2*V_pos(:,i) + V_pos(:,i+1)).^2;
            IS_pos_3 = 1/4 * (3*V_pos(:,i) - 4*V_pos(:,i+1)+ V_pos(:,i+2)).^2 + 13/12 * (V_pos(:,i) - 2*V_pos(:,i+1) + V_pos(:,i+2)).^2;

            alpha_pos_1 = 0; 
            alpha_pos_2 = 6/10./ (ones(3,1) * epsilon + IS_pos_2).^2; 
            alpha_pos_3 = 3/10./ (ones(3,1) * epsilon + IS_pos_3).^2; 

            omega_pos_1 = alpha_pos_1 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);
            omega_pos_2 = alpha_pos_2 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);
            omega_pos_3 = alpha_pos_3 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);

            V_pos_1 = 0;
            V_pos_2 = -1/6 * V_pos(:,i-1) + 5/6 * V_pos(:,i) + 1/3 * V_pos(:,i+1);
            V_pos_3 = 1/3 * V_pos(:,i) + 5/6 * V_pos(:,i+1) - 1/6 * V_pos(:,i+2);
            % j-1/2
            IS_neg_1 = 1/4 * (V_neg(:,i+2) - 4*V_neg(:,i+1) + 3*V_neg(:,i)).^2 + 13/12 * (V_neg(:,i+2) - 2*V_neg(:,i+1) + V_neg(:,i)).^2;
            IS_neg_2 = 1/4 * (V_neg(:,i+1) - V_neg(:,i-1)).^2 + 13/12 * (V_neg(:,i+1) - 2*V_neg(:,i) + V_neg(:,i-1)).^2;

            alpha_neg_1 = 1/10./ (ones(3,1) * epsilon + IS_neg_1).^2; 
            alpha_neg_2 = 6/10./ (ones(3,1) * epsilon + IS_neg_2).^2; 
            alpha_neg_3 = 0; 

            omega_neg_1 = alpha_neg_1 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);
            omega_neg_2 = alpha_neg_2 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);
            omega_neg_3 = alpha_neg_3 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);

            V_neg_1 = 1/3 * V_neg(:,i+2) - 7/6 * V_neg(:,i+1) + 11/6 * V_neg(:,i);
            V_neg_2 = -1/6 * V_neg(:,i+1) + 5/6 * V_neg(:,i) + 1/3 * V_neg(:,i-1);
            V_neg_3 = 0;
        elseif i == N+2
            for k = i-2:i+1
                V_pos(:,k) = lambda_pos * S * U(:,k) ;
                V_neg(:,k) = lambda_neg * S * U(:,k) ;
            end
            % WENO
            % j+1/2
            IS_pos_1 = 1/4 * (V_pos(:,i-2) - 4*V_pos(:,i-1)+ 3*V_pos(:,i)).^2 + 13/12 * (V_pos(:,i-2) - 2*V_pos(:,i-1) + V_pos(:,i)).^2;
            IS_pos_2 = 1/4 * (V_pos(:,i-1) - V_pos(:,i+1)).^2 + 13/12 * (V_pos(:,i-1) - 2*V_pos(:,i) + V_pos(:,i+1)).^2;

            alpha_pos_1 = 1/10./ (ones(3,1) * epsilon + IS_pos_1).^2; 
            alpha_pos_2 = 6/10./ (ones(3,1) * epsilon + IS_pos_2).^2; 
            alpha_pos_3 = 0; 

            omega_pos_1 = alpha_pos_1 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);
            omega_pos_2 = alpha_pos_2 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);
            omega_pos_3 = alpha_pos_3 ./ (alpha_pos_1 + alpha_pos_2 + alpha_pos_3);

            V_pos_1 = 1/3 * V_pos(:,i-2) - 7/6 * V_pos(:,i-1) + 11/6 * V_pos(:,i);
            V_pos_2 = -1/6 * V_pos(:,i-1) + 5/6 * V_pos(:,i) + 1/3 * V_pos(:,i+1);
            V_pos_3 = 0;
            % j-1/2
            IS_neg_2 = 1/4 * (V_neg(:,i+1) - V_neg(:,i-1)).^2 + 13/12 * (V_neg(:,i+1) - 2*V_neg(:,i) + V_neg(:,i-1)).^2;
            IS_neg_3 = 1/4 * (3*V_neg(:,i) - 4*V_neg(:,i-1)+ V_neg(:,i-2)).^2 + 13/12 * (V_neg(:,i) - 2*V_neg(:,i-1) + V_neg(:,i-2)).^2;

            alpha_neg_1 = 0; 
            alpha_neg_2 = 6/10./ (ones(3,1) * epsilon + IS_neg_2).^2; 
            alpha_neg_3 = 3/10./ (ones(3,1) * epsilon + IS_neg_3).^2; 

            omega_neg_1 = alpha_neg_1 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);
            omega_neg_2 = alpha_neg_2 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);
            omega_neg_3 = alpha_neg_3 ./ (alpha_neg_1 + alpha_neg_2 + alpha_neg_3);

            V_neg_1 = 0;
            V_neg_2 = -1/6 * V_neg(:,i+1) + 5/6 * V_neg(:,i) + 1/3 * V_neg(:,i-1);
            V_neg_3 = 1/3 * V_neg(:,i) + 5/6 * V_neg(:,i-1) - 1/6 * V_neg(:,i-2);
        end    
        V_pos_c(:,i) = omega_pos_1 .* V_pos_1 + omega_pos_2 .* V_pos_2 + omega_pos_3 .* V_pos_3;    % j+1/2
        V_neg_c(:,i) = omega_neg_1 .* V_neg_1 + omega_neg_2 .* V_neg_2 + omega_neg_3 .* V_neg_3;    % j-1/2
    end
 
    for i = 1:N+1 
       S_inv = [1                           1                      1;
                u_bar(i)-a_bar(i)            u_bar(i)               u_bar(i)+a_bar(i);
                H_bar(i)-u_bar(i)*a_bar(i)   1/2*u_bar(i)^2         H_bar(i)+u_bar(i)*a_bar(i)] ;
       F_c(:,i) = S_inv * (V_pos_c(:,i) + V_neg_c(:,i+1));
    end
    F_c(:,N+2) = F_c(:,N+1);
    for i = 2:N+2      
        L(:,i) = - 1/dx * (F_c(:,i) - F_c(:,i-1));
    end
    L(:,1) = 0;      %transmissive left boundary
    L(:,N+3) = 0;  %transmissive right boundary
end
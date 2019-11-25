function [rho, u, p, e] = roe_solver(dx, rho0, u0, p0, tEnd, entropy_fix)
%input initial condition & total time
global gamma CFL N epsilon
t = 0;
rho = rho0;
u = u0;
p = p0;
E = p ./ ((gamma-1) * rho) + u.^2/2;
U = [rho; rho.*u; rho.*E];   % U values at i
entropy_fix = input('please input entropy fix type: \n  {0} no entropy fix; {1} entropy fix depend on epsilon; {2} Harten-Hyman entropy fix. \n ');
while t <= tEnd
    % update time step using CFL condition
    a = sqrt(gamma * p ./ rho);
    S_max = max(max(abs(u) + a));
    dt = dx * CFL / S_max;
    % time advance
    U_sub = U;
    L = Calculate_L_Roe(U, entropy_fix, dx);
    U = U_sub + dt * L;
    % compute flow properties
    rho = U(1,:);
    u = U(2,:)./rho;
    E = U(3,:)./rho;
    p = (gamma-1).*(E-u.^2./2).*rho;
    e = p./(gamma-1)./rho;
    t = t+dt;
end
end
function L = Calculate_L_Roe(U, entropy_fix, dx)
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
    
    if entropy_fix ==1       % entropy fix
        for i = 1:N+2
            if abs(lambda_bar_1(i)) <= epsilon
                lambda_bar_1(i) = (lambda_bar_1(i)^2 + epsilon^2) / 2 / epsilon;
            end
            if abs(lambda_bar_2(i)) <= epsilon
                lambda_bar_2(i) = (lambda_bar_2(i)^2 + epsilon^2) / 2 / epsilon;
            end
            if abs(lambda_bar_3(i)) <= epsilon
                lambda_bar_3(i) = (lambda_bar_3(i)^2 + epsilon^2) / 2 / epsilon;
            end
        end 
    elseif entropy_fix == 2     % Harten-Hyman entropy fix
        lambda_1 = u - a;
        lambda_2 = u;
        lambda_3 = u + a;
        for i = 1:N+2
            epsilon_1 = max([0, (lambda_bar_1(i)-lambda_1(i)), (lambda_1(i+1)-lambda_bar_1(i))]);
            if abs(lambda_bar_1(i)) < epsilon_1
                lambda_bar_1(i) = epsilon_1;
            else
                lambda_bar_1(i) = abs(lambda_bar_1(i));
            end
            epsilon_2 = max([0, (lambda_bar_2(i)-lambda_2(i)), (lambda_2(i+1)-lambda_bar_2(i))]);
            if abs(lambda_bar_2(i)) < epsilon_2
                lambda_bar_2(i) = epsilon_2;
            else
                lambda_bar_2(i) = abs(lambda_bar_2(i));
            end
            epsilon_3 = max([0, (lambda_bar_3(i)-lambda_3(i)), (lambda_3(i+1)-lambda_bar_3(i))]);
            if abs(lambda_bar_3(i)) < epsilon_3
                lambda_bar_3(i) = epsilon_3;
            else
                lambda_bar_3(i) = abs(lambda_bar_3(i));
            end
        end 
    end
    lambda_bar_1=repmat(lambda_bar_1,3,1);
    lambda_bar_2=repmat(lambda_bar_2,3,1);
    lambda_bar_3=repmat(lambda_bar_3,3,1);
    %averaged right eigenvectors
    K_bar_1 = [ones(1,N+2); u_bar-a_bar; H_bar-u_bar.*a_bar];
    K_bar_2 = [ones(1,N+2); u_bar; u_bar.^2/2];
    K_bar_3 = [ones(1,N+2); u_bar+a_bar; H_bar+u_bar.*a_bar];
    %wave strengths
    Delta_U = U(:,2:N+3)-U(:,1:N+2);
    alpha_bar_2 = (gamma-1) ./ (a_bar.^2) .* (H_bar-u_bar.^2) .* Delta_U(1,:)...
                + (gamma-1) .* u_bar ./ a_bar.^2 .* Delta_U(2,:)...
                - (gamma-1) ./ a_bar.^2 .* Delta_U(3,:);
    alpha_bar_1 = ones(1,N+2) ./ (2.*a_bar) .* (Delta_U(1,:) .* (u_bar+a_bar) - Delta_U(2,:) - a_bar.*alpha_bar_2);
    alpha_bar_3 = Delta_U(1,:) - (alpha_bar_1 + alpha_bar_2);
    alpha_bar_1 = repmat(alpha_bar_1,3,1);
    alpha_bar_2 = repmat(alpha_bar_2,3,1);
    alpha_bar_3 = repmat(alpha_bar_3,3,1);
    
    F = [rho.*u; rho.*u.^2+p; rho.*u.*H]; % F(U_j^n)
    % numerical Roe flux
    RF = 0.5 * (F(:,1:N+2)+F(:,2:N+3))...
        - (abs(lambda_bar_1).*alpha_bar_1.*K_bar_1...
        + abs(lambda_bar_2).*alpha_bar_2.*K_bar_2...
        + abs(lambda_bar_3).*alpha_bar_3.*K_bar_3) * 0.5; %flux acorss cells

    % calculate L
    for i = 2:N+2
        L(:,i) = 1/dx * (RF(:,i-1) - RF(:,i));
    end
    L(:,1) = 0;      %transmissive left boundary
    L(:,N+3) = 0;  %transmissive right boundary
end
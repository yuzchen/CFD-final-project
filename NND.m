function [rho, u, p, e] = NND(dx, rho0, u0, p0, tEnd)
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
    % 1st step
    U_sub = U;
    L = Calculate_L_NND(U, dx);
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
function  L = Calculate_L_NND(U, dx)
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
    % calculate flux
    F_c = zeros(1,N+2);
    for k = 1:3
        for i = 2:N+1
            F_c(k,i) = F_pos(k,i) + F_neg(k,i+1) + 1/2 * minmod([F_pos(k,i)-F_pos(k,i-1), F_pos(k,i+1)-F_pos(k,i)]) - 1/2 * minmod([F_neg(k,i+1)-F_neg(k,i), F_neg(k,i+2)-F_neg(k,i+1)]);
        end
    end
    F_c(:,1) = F_c(:,2); F_c(:,N+2) = F_c(:,N+1);
    for i = 2:N+2
        L(:,i) = 1/dx * (F_c(:,i-1) - F_c(:,i));
    end
    L(:,1) = 0;      %transmissive left boundary
    L(:,N+3) = 0;  %transmissive right boundary
end
function [rho, u, p, e] = MacCormack(dx, rho0, u0, p0, tEnd)
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
    [U_c, L] = Calculate_L_MC(U, dx, dt);
    U_c(:,N+3) = U_c(:,N+2);
    U = 1/2 * (U_c + U) + dt/2 * L;   
    % compute flow properties
    rho = U(1,:);
    u = U(2,:)./rho;
    E = U(3,:)./rho;
    p = (gamma-1).*(E-u.^2./2).*rho;
    e = p./(gamma-1)./rho;
    t = t+dt;
end
end
function [U_c, L] = Calculate_L_MC(U, dx, dt)
    global gamma CFL N epsilon
    rho = U(1,:);
    u = U(2,:)./rho;
    E = U(3,:)./rho;
    p = (gamma-1).*(E-u.^2./2).*rho;
    H = (p * gamma) ./ ((gamma-1) * rho) + u.^2/2;
    a = sqrt(gamma * p ./ rho);
    F = [rho.*u; rho.*u.^2+p; rho.*u.*H]; % F(U_j^n)
    U_c = U(:,1:N+2) - dt/dx * (F(:,2:N+3) - F(:,1:N+2));  % U_{i}^{n+1/2}
    rho_c = U_c(1,:);
    u_c = U_c(2,:)./rho_c;
    E_c = U_c(3,:)./rho_c;
    p_c = (gamma-1).*(E_c-u_c.^2./2).*rho_c;
    H_c = (p_c * gamma) ./ ((gamma-1) * rho_c) + u_c.^2/2;
    F = [rho_c.*u_c; rho_c.*u_c.^2+p_c; rho_c.*u_c.*H_c]; % F(U_{i}^{n+1/2})
    % calculate L
    for i = 2:N+2
        L(:,i) = 1/dx * (F(:,i-1) - F(:,i));
    end
    L(:,1) = 0;      %transmissive left boundary
    L(:,N+3) = 0;  %transmissive right boundary
end
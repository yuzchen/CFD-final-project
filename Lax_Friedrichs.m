function [rho, u, p, e] = Lax_Friedrichs(dx, rho0, u0, p0, tEnd)
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
    L = Calculate_L_LF(U, dx, dt);
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
function L = Calculate_L_LF(U, dx, dt)
    global gamma CFL N epsilon
    rho = U(1,:);
    u = U(2,:)./rho;
    E = U(3,:)./rho;
    p = (gamma-1).*(E-u.^2./2).*rho;
    H = (p * gamma) ./ ((gamma-1) * rho) + u.^2/2;
    a = sqrt(gamma * p ./ rho);
    F = [rho.*u; rho.*u.^2+p; rho.*u.*H]; % F(U_j^n)
    % numerical flux
    LF = 0.5 * (F(:,1:N+2)+F(:,2:N+3) - dx / dt * (U(:,2:N+3) - U(:,1:N+2)));
    % calculate L
    for i = 2:N+2
        L(:,i) = 1/dx * (LF(:,i-1) - LF(:,i));
    end
    L(:,1) = 0;      %transmissive left boundary
    L(:,N+3) = 0;  %transmissive right boundary
end
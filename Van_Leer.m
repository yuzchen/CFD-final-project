function [rho, u, p, e] = Van_Leer(dx, rho0, u0, p0, tEnd, time_advance_type)
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
    if time_advance_type == 1
%         dt = dx * CFL / S_max;
        dt = 0.0001;
    elseif time_advance_type == 2
        dt = 0.0001;
    end
    % time advance
    U_sub = U;
    L = Calculate_L_VL(U, dx, time_advance_type);
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
function  L = Calculate_L_VL(U, dx, time_advance_type)
    global gamma CFL N epsilon
    rho = U(1,:);
    u = U(2,:)./rho;
    E = U(3,:)./rho;
    p = (gamma-1).*(E-u.^2./2).*rho;
    H = (p * gamma) ./ ((gamma-1) * rho) + u.^2/2;
    a = sqrt(gamma * p ./ rho);
    Ma = u ./ a;
    F = [rho.*u; rho.*u.^2+p; rho.*u.*H]; 
    F_1_pos = rho .* a .* ((Ma + 1)/2).^2;
    F_1_neg = - rho .* a .* ((Ma - 1)/2).^2;
    for i = 1:N+3
        if Ma(i) >= 1
            F_pos(:,i) = F(:,i);
            F_neg = zeros(3,N+3);
        elseif  Ma(i) >= -1
            F_pos(:,i) = [F_1_pos(i); 
                          F_1_pos(i)/gamma * ((gamma-1)*u(i) + 2*a(i));
                          F_1_pos(i)/2/(gamma^2-1) * ((gamma-1)*u(i) + 2*a(i))^2];
            F_neg(:,i) = [F_1_neg(i); 
                          F_1_neg(i)/gamma * ((gamma-1)*u(i) - 2*a(i));
                          F_1_neg(i)/2/(gamma^2-1) * ((gamma-1)*u(i) - 2*a(i))^2];          
        else
            F_neg(:,i) = F(:,i);
            F_pos = zeros(3,N+3);
        end
    end
    % apply upwind scheme
    if time_advance_type == 1   % 1st order
        for i = 2:N+2
            L(:,i) = - 1/dx * (F_pos(:,i) - F_pos(:,i-1) + F_neg(:,i+1) - F_neg(:,i));
        end
        L(:,1) = -1/dx*(F_neg(:,2) - F_neg(:,1));      
        L(:,N+3) = -1/dx*(F_pos(:,N+3) - F_pos(:,N+2));  
    elseif time_advance_type ==2    % 2nd order
        for i = 3:N+1
            L(:,i) = - 1/2/dx * (3 * F_pos(:,i) - 4 * F_pos(:,i-1) + F_pos(:,i-2) - F_neg(:,i+2) + 4*F_neg(:,i+1) - 3*F_neg(:,i));
        end
        L(:,2) = - 1/dx * (F_pos(:,2) - F_pos(:,1) + F_neg(:,3) - F_neg(:,2));
        L(:,N+2) = - 1/dx * (F_pos(:,N+2) - F_pos(:,N+1) + F_neg(:,N+3) - F_neg(:,N+2));
        L(:,1) = -1/dx*(F_neg(:,2) - F_neg(:,1));      
        L(:,N+3) = -1/dx*(F_pos(:,N+3) - F_pos(:,N+2));  
    end
end
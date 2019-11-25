function [rho, u, p, e] = LBM_KT(dx, rho0, u0, p0, tEnd, time_advance_type)   % KT model
    global gamma CFL N epsilon
    t = 0; 
    rho = rho0;
    u = u0;
    p = p0;
    E = p ./ ((gamma-1) * rho) + u.^2/2;
    U = [rho; rho.*u; rho.*E];   % U values at j
    knudsen = 0.0001; dt = knudsen / 4;
    v1 = 1;     v2 = 3;    eta0 = 2;  % free parameters
    b = 2 / (gamma - 1);
    eta_i = zeros(5,N+3);
    eta_i(1,:) = eta0 * ones(1,N+3);
    c_i = [zeros(1,N+3);
            v1 * ones(1,N+3);
           -v1 * ones(1,N+3);
            v2 * ones(1,N+3);
           -v2 * ones(1,N+3)];
    f = calculate_feq(v1, v2, eta0, b, c_i, U);  % initial value
    L = zeros(5,N+3);
    while t <= tEnd
        feq = calculate_feq(v1, v2, eta0, b, c_i, U);    % local equilibrium velocity distribution fucntion
        % solve DVBE 1rd in time and 1rd upwind in space
        f_pos = 1/2 .* (c_i + abs(c_i)) .* f;
        f_neg = 1/2 .* (c_i - abs(c_i)) .* f;
        if time_advance_type == 1    % 1st order
            for i = 2:N+2
                L(:,i) = - 1/dx * (f_pos(:,i) - f_pos(:,i-1) + f_neg(:,i+1) - f_neg(:,i));
            end
            L(:,1) = - 1/dx * (f_neg(:,2) - f_neg(:,1));      
            L(:,N+3) = - 1/dx * (f_pos(:,N+3) - f_pos(:,N+2));  
         elseif time_advance_type == 2    % 2nd order
            for i = 3:N+1
                L(:,i) = - 1/2/dx * (3 * f_pos(:,i) - 4 * f_pos(:,i-1) + f_pos(:,i-2) - f_neg(:,i+2) + 4*f_neg(:,i+1) - 3*f_neg(:,i));
            end
            L(:,2) = - 1/dx * (f_pos(:,2) - f_pos(:,1) + f_neg(:,3) - f_neg(:,2));
            L(:,N+2) = - 1/dx * (f_pos(:,N+2) - f_pos(:,N+1) + f_neg(:,N+3) - f_neg(:,N+2));
            L(:,1) = -1/dx*(f_neg(:,2) - f_neg(:,1));      
            L(:,N+3) = -1/dx*(f_pos(:,N+3) - f_pos(:,N+2));  
        end
        f = f + dt * L + (feq - f) / knudsen * dt;
        U = macroscopic_variables(c_i, eta_i, f);
        t = t + dt;            
    end
    rho = U(1,:);
    u = U(2,:)./rho;
    E = U(3,:)./rho;
    p = (gamma-1).*(E-u.^2./2).*rho;
    e = p./(gamma-1)./rho;
end

function f = calculate_feq(v1, v2, eta0, b, c_i, U)
    global gamma N
    A_i = zeros(5,N+3); B_i = zeros(5,N+3);
    rho = U(1,:);
    u = U(2,:) ./ U(1,:);
    T = (2 * U(3,:) - rho .* u.^2) ./ rho / b;
    A_i(1,:) = (b-1) / eta0^2 .* T;
    A_i(2:3,:) = repmat(1 / (2*(v1^2 - v2^2)) * ( - v2^2 * ones(1,N+3) + ((b-1)*v2^2/eta0^2 + 1) .* T + u.^2), 2, 1);
    A_i(4:5,:) = repmat(1 / (2*(v2^2 - v1^2)) * ( - v1^2 * ones(1,N+3) + ((b-1)*v1^2/eta0^2 + 1) .* T + u.^2), 2, 1);
    B_i(1,:) = zeros(1,N+3);
    B_i(2:3,:) = repmat((-v2^2 * ones(1,N+3) + (b+2).*T + u.^2) / (2 * v1^2 * (v1^2 - v2^2)), 2, 1);
    B_i(4:5,:) = repmat((-v1^2 * ones(1,N+3) + (b+2).*T + u.^2) / (2 * v2^2 * (v2^2 - v1^2)), 2, 1);
    f = rho .* (A_i + B_i .* u .* c_i);
end

function U = macroscopic_variables(c_i, eta_i, f)
    U(1,:) = sum(f,1);
    U(2,:) = sum(f .* c_i, 1);
    U(3,:) = sum(f .* (c_i.^2 + eta_i.^2)) / 2;
end

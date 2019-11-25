function [rho, u, p, e] = LBM_Qu(dx, rho0, u0, p0, tEnd)   % Qu's model D1Q5L2
    global gamma CFL N epsilon
    t = 0; 
    rho = rho0;
    u = u0;
    p = p0;
    E = p ./ ((gamma-1) * rho) + u.^2/2;
    U = [rho; rho.*u; rho.*E];   % U values at j
    knudsen = 0.0001; dt = knudsen / 4;
    d1 = 1;    d2 = 2;   lambda2 = 1;% parameters
    b = 2 / (gamma - 1);
    e_i = [zeros(1,N+3);
            d1 * ones(1,N+3);
           -d1 * ones(1,N+3);
            d2 * ones(1,N+3);
           -d2 * ones(1,N+3)];
    [fi1, fi2] = calculate_feq(d1, d2, lambda2, b, U);
    L_1 = zeros(5,N+3); L_2 = zeros(5,N+3);
    while t <= tEnd
        [feq1, feq2] = calculate_feq(d1, d2, lambda2, b, U);
        % solve DVBE 1rd in time and 1rd upwind in space
        f_1_pos = 1/2 .* (e_i + abs(e_i)) .* fi1;
        f_1_neg = 1/2 .* (e_i - abs(e_i)) .* fi1;
        f_2_pos = 1/2 .* (e_i + abs(e_i)) .* fi2;
        f_2_neg = 1/2 .* (e_i - abs(e_i)) .* fi2;
        for i = 2:N+2
            L_1(:,i) = - 1/dx * (f_1_pos(:,i) - f_1_pos(:,i-1) + f_1_neg(:,i+1) - f_1_neg(:,i));
            L_2(:,i) = - 1/dx * (f_2_pos(:,i) - f_2_pos(:,i-1) + f_2_neg(:,i+1) - f_2_neg(:,i));
        end
        L_1(:,1) = - 1/dx * (f_1_neg(:,2) - f_1_neg(:,1));      
        L_2(:,1) = - 1/dx * (f_2_neg(:,2) - f_2_neg(:,1));      
        L_1(:,N+3) = - 1/dx * (f_1_pos(:,N+3) - f_1_pos(:,N+2));  
        L_2(:,N+3) = - 1/dx * (f_2_pos(:,N+3) - f_2_pos(:,N+2));  
        fi1 = fi1 + dt * L_1 + (feq1 - fi1) / knudsen * dt;
        fi2 = fi2 + dt * L_2 + (feq2 - fi2) / knudsen * dt;
        U = macroscopic_variables(lambda2, e_i, fi1, fi2);
        t = t + dt;            
    end
    rho = U(1,:);
    u = U(2,:)./rho;
    E = U(3,:)./rho;
    p = (gamma-1).*(E-u.^2./2).*rho;
    e = p./(gamma-1)./rho;
end

function [fi1, fi2] = calculate_feq(d1, d2, lambda2, b, U)
    global gamma N
    rho_i = zeros(5,N+3); 
    rho = U(1,:);
    u = U(2,:) ./ U(1,:);
    E = U(3,:)./rho;
    p = (gamma-1).*(E-u.^2./2).*rho;
    e = p./(gamma-1)./rho;
    c = sqrt((gamma-1).*e);
    e_p = (3-gamma)/2 .* e;
    rho_i(1,:) = rho .* (d1^2 * d2^2 * ones(1,N+3) - d1^2 * u.^2 - d1^2 * c.^2 -  d2^2 * u.^2 - d2^2 * c.^2 + u.^4 + 6*u.^2 .* c.^2 + c.^4) / (d1^2 * d2^2);
    rho_i(2,:) = rho .* (3 * d1 .* u .* c.^2 - d1 * d2^2 * u - d2^2 * u.^2 - d2^2 * c.^2 + 6*u.^2 .* c.^2 + c.^4 + u.^4 + d1 * u.^3) / (2 * d1^2 * (d1^2 - d2^2));
    rho_i(3,:) = rho .* (-3 * d1 .* u .* c.^2 + d1 * d2^2 * u - d2^2 * u.^2 - d2^2 * c.^2 + 6*u.^2 .* c.^2 + c.^4 + u.^4 - d1 * u.^3) / (2 * d1^2 * (d1^2 - d2^2));
    rho_i(4,:) = - rho .* (-d2*d1^2*u + 3*d2*u.*c.^2 - d1^2*u.^2 - d1^2.*c.^2 + 6*u.^2.*c.^2 + c.^4 + u.^4 + d2*u.^3) / (2 * d2^2 * (d1^2 - d2^2));
    rho_i(5,:) = - rho .* (d2*d1^2*u - 3*d2*u.*c.^2 - d1^2*u.^2 - d1^2.*c.^2 + 6*u.^2.*c.^2 + c.^4 + u.^4 - d2*u.^3) / (2 * d2^2 * (d1^2 - d2^2));
    fi1 = rho_i .* (lambda2 * ones(5,N+3) - repmat(e_p,5,1)) / lambda2;
    fi2 = rho_i .* repmat(e_p,5,1) / lambda2;
end

function U = macroscopic_variables(lambda2, e_i, fi1, fi2)
    global N
    U(1,:) = sum(fi1,1) + sum(fi2,1);
    U(2,:) = sum(fi1 .* e_i + fi2 .* e_i, 1);
    U(3,:) = sum(((fi1+fi2) .* e_i.^2 + 2 * fi2 .* lambda2), 1) / 2;
end

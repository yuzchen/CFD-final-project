function [rho_ana, u_ana, p_ana, e_ana] = analytical_solution(xn, x0, rho0, u0, p0, tEnd)
global gamma N
xn = xn - x0;
rho1 = rho0(1);
rho2 = rho0(end);
u1 = u0(1);
u2 = u0(end);
p1 = p0(1);
p2 = p0(end);
c1 = sqrt(gamma*p1 / rho1);
c2 = sqrt(gamma*p2 / rho2);
rho_ana = zeros(size(xn));
u_ana = zeros(size(xn));
p_ana = zeros(size(xn));
if p0(end) >= p0(1)
    F_p2 = (p2 - p1) / (rho1 * c1 * sqrt((gamma + 1) / 2 / gamma * p2 / p1 + (gamma - 1) / 2 /gamma));
    F_p1 = - 2 * c2 / (gamma -1) * (1 - (p1 / p2)^((gamma-1)/2/gamma));
    F_0  = - 2 * c1 / (gamma - 1) - 2 * c2 / (gamma - 1);
    if (u1 - u2) >= F_p2    %left and right shock wave
        disp(['solution type = 1']);           
        f = @(x)((x - p1) / rho1 / c1 / sqrt((gamma+1) / 2 / gamma * x / p1 + (gamma-1) / 2 / gamma)...
               + (x - p2) / rho2 / c2 / sqrt((gamma+1) / 2 / gamma * x / p2 + (gamma-1) / 2 / gamma) - (u1 - u2));
        p_bar = fsolve(f, 0);
        A1 = rho1 * c1 * sqrt((gamma+1) / gamma / 2 * p_bar / p1 + (gamma - 1) / 2 /gamma);
        A2 = rho2 * c2 * sqrt((gamma+1) / gamma / 2 * p_bar / p2 + (gamma - 1) / 2 /gamma);
        U = 0.5 * (u1 + u2 + (p_bar - p2) / A2 - (p_bar - p1) / A1);
        Z1 = u1 - A1 / rho1;
        R1 = rho1 * A1 / (A1 - rho1 * (u1 - U));       
        Z2 = u2 + A2 / rho2;
        R2 = rho2 * A2 / (A2 + rho2 * (u2 - U));   
        for i = 1:N+2
            if (xn(i) / tEnd < Z1)
                rho_ana(i) = rho1;
                u_ana(i) = u1;
                p_ana(i) = p1;
            elseif (xn(i) / tEnd < U)
                rho_ana(i) = R1;
                u_ana = U;
                p_ana = p_bar;
            elseif (xn(i) / tEnd < Z2)
                rho_ana(i) = R2;
                u_ana(i) = U;
                p_ana(i) = p_bar;
            else 
                rho_ana(i) = rho2;
                u_ana(i) = u2;
                p_ana(i) = p2;
            end
        end             
    elseif ((u1 - u2) < F_p2) && ((u1 - u2) >= F_p1)   % left shock wave right rarefaction wave
        disp(['solution type = 3']);  
        f = @(x)(x - p1) / rho1 / c1 / sqrt((gamma+1) / 2 / gamma * x / p1 + (gamma-1) / 2 / gamma)...
               + 2 * c2 / (gamma - 1) * ((x / p2)^((gamma-1)/2/gamma) - 1) - (u1 - u2);
        p_bar = fsolve(f, 0);
        A1 = rho1 * c1 * sqrt((gamma+1) / gamma / 2 * p_bar / p1 + (gamma - 1) / 2 /gamma);
        A2 = rho2 * c2 * sqrt((gamma+1) / gamma / 2 * p_bar / p2 + (gamma - 1) / 2 /gamma);
        U = 0.5 * (u1 + u2 + (p_bar - p2) / A2 - (p_bar - p1) / A1);
        Z1 = u1 - A1 / rho1;
        R1 = rho1 * A1 / (A1 - rho1 * (u1 - U));   
        c2_star = c2 - (gamma-1) / 2 * (u2 - U);
        Z2 = u2 + c2;
        Z2_star = U + c2_star;
        R2 = gamma * p_bar / c2_star^2;
        for i = 1:N+2
            if (xn(i) / tEnd < Z1)
                rho_ana(i) = rho1;
                u_ana(i) = u1;
                p_ana(i) = p1;
            elseif (xn(i) / tEnd < U)
                rho_ana(i) = R1;
                u_ana(i) = U;
                p_ana(i) = p_bar;
            elseif (xn(i) / tEnd < Z2_star)
                rho_ana(i) = R2;
                u_ana(i) = U;
                p_ana(i) = p_bar;
            elseif (xn(i) / tEnd < Z2)
                c = (gamma - 1) / (gamma + 1) * (xn(i) / tEnd - u2) + 2 / (gamma - 1) * c2;
                u_ana(i) = xn(i) / tEnd - c;
                p_ana(i) = p2 * (c / c2)^(2*gamma / (gamma-1));
                rho_ana(i) = gamma * p_ana(i) / c^2;
            else
                rho_ana(i) = rho2;
                u_ana(i) = u2;
                p_ana(i) = p2;
            end
        end  
    elseif ((u1 - u2) < F_p1) && ((u1 - u2) >= F_0)    % left and right rarefaction wave
        disp(['solution type = 4']);  
        f = @(x)2 * c1 / (gamma - 1) * ((x / p1)^((gamma-1)/2/gamma) - 1)...
              + 2 * c2 / (gamma - 1) * ((x / p2)^((gamma-1)/2/gamma) - 1) - (u1 - u2);
        p_bar = fsolve(f, 0);
        A1 = rho1 * c1 * sqrt((gamma+1) / gamma / 2 * p_bar / p1 + (gamma - 1) / 2 /gamma);
        A2 = rho2 * c2 * sqrt((gamma+1) / gamma / 2 * p_bar / p2 + (gamma - 1) / 2 /gamma);
        U = 0.5 * (u1 + u2 + (p_bar - p2) / A2 - (p_bar - p1) / A1);
        c1_star = c1 + (gamma-1) / 2 * (u1 - U);
        Z1 = u1 - c1;
        Z1_star = U - c1_star;
        R1 = gamma * p_bar / c1_star^2;
        c2_star = c2 - (gamma-1) / 2 * (u2 - U);
        Z2 = u2 + c2;
        Z2_star = U + c2_star;
        R2 = gamma * p_bar / c2_star^2;
        for i = 1:N+2
            if (xn(i) / tEnd < Z1)
                rho_ana(i) = rho1;
                u_ana(i) = u1;
                p_ana(i) = p1;
            elseif (xn(i) / tEnd < Z1_star)
                c = (gamma - 1) / (gamma + 1) * (u1 - xn(i) / tEnd) + 2 / (gamma - 1) * c1;
                u_ana(i) = xn(i) / tEnd + c;
                p_ana(i) = p1 * (c / c1)^(2*gamma / (gamma-1));
                rho_ana(i) = gamma * p_ana(i) / c^2;
            elseif (xn(i) / tEnd < U)
                rho_ana(i) = R1;
                u_ana(i) = U;
                p_ana(i) = p_bar;
            elseif (xn(i) / tEnd < Z2_star)
                rho_ana(i) = R2;
                u_ana(i) = U;
                p_ana(i) = p_bar;
            elseif (xn(i) / tEnd < Z2)
                c = (gamma - 1) / (gamma + 1) * (xn(i) / tEnd - u2) + 2 / (gamma - 1) * c2;
                u_ana(i) = xn(i) / tEnd - c;
                p_ana(i) = p2 * (c / c2)^(2*gamma / (gamma-1));
                rho_ana(i) = gamma * p_ana(i) / c^2;
            else
                rho_ana(i) = rho2;
                u_ana(i) = u2;
                p_ana(i) = p2;
            end
        end  
    elseif (u1 - u2) < F_0                             % left and right rarefaction and vacuum
        disp(['solution type = 5, RCVCR']);  
        p_bar = 0;
        Z1 = u1 - c1;
        Z1_star = u1 + 2 / (gamma - 1) * c1;
        Z2 = u2 + c2;
        Z2_star = u2 - 2 / (gamma - 1) * c2;
        for i = 1:N+2
            if (xn(i) / tEnd < Z1)
                rho_ana(i) = rho1;
                u_ana(i) = u1;
                p_ana(i) = p1;
            elseif (xn(i) / tEnd < Z1_star)
                c = (gamma - 1) / (gamma + 1) * (u1 - xn(i) / tEnd) + 2 / (gamma - 1) * c1;
                u_ana(i) = xn(i) / tEnd + c;
                p_ana(i) = p1 * (c / c1)^(2*gamma / (gamma-1));
                rho_ana(i) = gamma * p_ana(i) / c^2;
            elseif (xn(i) / tEnd < Z2_star)
                rho_ana(i) = 0;
                u_ana(i) = 0;
                p_ana(i) = p_bar;
            elseif (xn(i) / tEnd < Z2)
                c = (gamma - 1) / (gamma + 1) * (xn(i) / tEnd - u2) + 2 / (gamma - 1) * c2;
                u_ana(i) = xn(i) / tEnd - c;
                p_ana(i) = p2 * (c / c2)^(2*gamma / (gamma-1));
                rho_ana(i) = gamma * p_ana(i) / c^2;
            else
                rho_ana(i) = rho2;
                u_ana(i) = u2;
                p_ana(i) = p2;
            end
        end  
    end
else
    disp(['solution type = 2, RCS']);  
    options = optimset('MaxFunEvals',1e5);
    f = @(x)(x - p2) / rho2 / c2 / sqrt((gamma+1) / 2 / gamma * x / p2 + (gamma-1) / 2 / gamma)...
           + 2 * c1 / (gamma - 1) * ((x / p1)^((gamma-1)/2/gamma) - 1) - (u1 - u2);
    [p_bar, fval, exitflag] = fsolve(f,0,options);
    A2 = rho2 * c2 * sqrt((gamma+1) / gamma / 2 * p_bar / p2 + (gamma - 1) / 2 / gamma);
    U = 0.5 * (u1 + u2 + (p_bar - p2) / A2 - 2 * c1 / (gamma-1) * ((p_bar/p1)^((gamma-1)/2/gamma) - 1));
    c1_star = c1 + (gamma-1) / 2 * (u1 - U);
    Z1 = u1 - c1;
    R1 = gamma * p_bar / c1_star^2;   
    Z1_star = U - c1_star;
    Z2 = u2 + A2 / rho2;
    R2 = rho2 * A2 / (A2 + rho2 * (u2 - U));
    for i = 1:N+2
        if (xn(i) / tEnd < Z1)
            rho_ana(i) = rho1;
            u_ana(i) = u1;
            p_ana(i) = p1;
        elseif (xn(i) / tEnd < Z1_star)
            c = (gamma-1) / (gamma+1) * (u1 - xn(i)/tEnd) + 2 / (gamma +1) * c1;
            u_ana(i) = xn(i)/tEnd + c;
            p_ana(i) = p1 * (c / c1) ^ (2 * gamma / (gamma -1));
            rho_ana(i) = gamma * p_ana(i) / c^2;
        elseif (xn(i) / tEnd < U)
            rho_ana(i) = R1;
            u_ana(i) = U;
            p_ana(i) = p_bar;
        elseif (xn(i) / tEnd < Z2)
            u_ana(i) = U;
            p_ana(i) = p_bar;
            rho_ana(i) = R2;
        else
            rho_ana(i) = rho2;
            u_ana(i) = u2;
            p_ana(i) = p2;
        end
    end  
end
e_ana = p_ana./(gamma-1)./rho_ana;
end
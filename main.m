clc; clear; close all;
%% define parameters
global gamma CFL N epsilon
gamma = 1.4;
CFL = 0.9;
N = 100;                        % number of cells; must be even
epsilon = 10^(-6);
initType = 1;                   % {1} Sod;  {3} Lax;
solverType = 14;                 % {1} Lax-Friedrichs;    {2} two-step Lax-Wendroff;   {3} MacCormack;    
                                % {4} Steger-Warming;    {5} Van Leer                 {6} L-F aplit;      {7} SCM;
                                % {8} Roe;               {9} Harten's TVD;            {10} NND; 
                                % {11} WENO-5-SW-split;  {12} WENO-5-SCM;
                                % {13} LBM-KT            {14} LBM-Qu;                
                                % {15} SPH, not working :(
upwind_order = 1;               % order of upwind scheme in case {4}{5}{6}{13}{14}
range = [-1 1];                 % spatial range from -1 to 1
%% geometry setup
dx = (range(2) - range(1)) / N;
xc = (range(1)-dx/2):dx:(range(2)+dx/2);                % centers location, including virtual points
xn = [range(1)-dx xc + dx/2];                           % nodes location
[rho0, u0, p0, tEnd, x0] = IC_setup(xn, initType);      % set up initial conditions at nodes
%% analytical solution
[rho_ana, u_ana, p_ana, e_ana] = analytical_solution(xn, x0, rho0, u0, p0, tEnd);
%% numerical solution
switch solverType
    case{1} % Lax-Friedrichs
        [rho_f, u_f, p_f, e_f] = Lax_Friedrichs(dx, rho0, u0, p0, tEnd);
    case{2} % Lax-Wendorff
        [rho_f, u_f, p_f, e_f] = Lax_Wendroff(dx, rho0, u0, p0, tEnd);
    case{3} % MacCormack
        [rho_f, u_f, p_f, e_f] = MacCormack(dx, rho0, u0, p0, tEnd);
    case{4} % Steger-Warming
        [rho_f, u_f, p_f, e_f] = Steger_Warming(dx, rho0, u0, p0, tEnd, upwind_order);
    case{5} % Van-Leer
        [rho_f, u_f, p_f, e_f] = Van_Leer(dx, rho0, u0, p0, tEnd, upwind_order);
    case{6} % L-F split
        [rho_f, u_f, p_f, e_f] = L_F_split(dx, rho0, u0, p0, tEnd, upwind_order);
    case{7} % SCM
        [rho_f, u_f, p_f, e_f] = SCM(dx, rho0, u0, p0, tEnd);
    case{8} % Roe_solver
        [rho_f, u_f, p_f, e_f] = roe_solver(dx, rho0, u0, p0, tEnd);
    case{9} % TVD
        [rho_f, u_f, p_f, e_f] = Harten_TVD(dx, rho0, u0, p0, tEnd);
    case{10} % NND
        [rho_f, u_f, p_f, e_f] = NND(dx, rho0, u0, p0, tEnd);
    case{11} % WENO-5-SW
        [rho_f, u_f, p_f, e_f] = WENO_SW(dx, rho0, u0, p0, tEnd);
    case{12} % WENO-5-SCM
        [rho_f, u_f, p_f, e_f] = WENO_SCM(dx, rho0, u0, p0, tEnd);
    case{13} % LBM-KT
        [rho_f, u_f, p_f, e_f] = LBM_KT(dx, rho0, u0, p0, tEnd, upwind_order);
    case{14} % LBM-Qu
        [rho_f, u_f, p_f, e_f] = LBM_Qu(dx, rho0, u0, p0, tEnd);
    case{15} % SPH
        SPH(xn, rho0, u0, p0, tEnd);
    otherwise 
        error('wrong solverType');
end
%% plot
figure(1);
subplot(2,2,1); plot(xn, rho_f,'-o','markersize',3); ylabel('Density'); xlabel('Position');
hold on; plot(xn,rho_ana,'r'); xlim(range);
subplot(2,2,2); plot(xn, u_f,'-o','markersize',3); xlim([0 1]); ylabel('Velocity'); xlabel('Position');
hold on; plot(xn,u_ana,'r'); xlim(range);
subplot(2,2,3); plot(xn, p_f,'-o','markersize',3); xlim([0 1]); ylabel('Pressure'); xlabel('Position');
hold on; plot(xn,p_ana,'r'); xlim(range);
subplot(2,2,4); plot(xn, e_f,'-o','markersize',3); xlim([0 1]); ylabel('Internal energy'); xlabel('Position');
hold on; plot(xn,e_ana,'r'); xlim(range);


% [rho_f1, u_f1, p_f1, e_f1] = WENO_SW(dx, rho0, u0, p0, tEnd);
% [rho_f2, u_f2, p_f2, e_f2] = Steger_Warming(dx, rho0, u0, p0, tEnd, 2);
 [rho_f1, u_f1, p_f1, e_f1] = Steger_Warming(dx, rho0, u0, p0, tEnd, 1);
 
 
hold on;subplot(2,2,1); plot(xn, rho_f1,'g-','markersize',3); ylabel('Density'); xlabel('Position');
hold on;  xlim(range);
subplot(2,2,2); plot(xn, u_f1,'g-','markersize',3); xlim([0 1]); ylabel('Velocity'); xlabel('Position');
hold on; xlim(range);
subplot(2,2,3); plot(xn, p_f1,'g-','markersize',3); xlim([0 1]); ylabel('Pressure'); xlabel('Position');
hold on; xlim(range);
subplot(2,2,4); plot(xn, e_f1,'g-','markersize',3); xlim([0 1]); ylabel('Internal energy'); xlabel('Position');
hold on; xlim(range);
% hold on;subplot(2,2,1); plot(xn, rho_f2,'m-','markersize',3); ylabel('Density'); xlabel('Position');
% hold on;  xlim(range);
% subplot(2,2,2); plot(xn, u_f2,'m-','markersize',3); xlim([0 1]); ylabel('Velocity'); xlabel('Position');
% hold on; xlim(range);
% subplot(2,2,3); plot(xn, p_f2,'m-','markersize',3); xlim([0 1]); ylabel('Pressure'); xlabel('Position');
% hold on; xlim(range);
% subplot(2,2,4); plot(xn, e_f2,'m-','markersize',3); xlim([0 1]); ylabel('Internal energy'); xlabel('Position');
% hold on; xlim(range);

% figure(2)
% plot(xn, rho_ana);hold on; plot(xn, u_ana); plot(xn, p_ana);xlim(range); grid on;




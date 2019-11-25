function [rho, v, p, e] = SPH(xn, rho_i, u0, p0, tEnd)
    global gamma
    m = 0.001875;
    N = 400;
    tEnd = 0.2;
    t = 0;
    dt = 0.001;
    for i = 1:320 % left
        rho0(i) = 1;   
        v(i) = 0;      
        p(i) = 1;      
        e(i) = 2.5; 
        x(i) = -0.6 + 0.6/320*i;
        h0(i) = 5*m/rho0(i)/4;  % initial smooth length
    end
    for i = 321:400   % right
        rho0(i) = 0.25;   
        v(i) = 0;      
        p(i) = 0.1795;     
        e(i) = 1.795; 
        x(i) = 0.6/80*(i-320);
        h0(i) = 5*m/rho0(i)/4;
    end
    rho = rho0;
    h = h0;
    e_m = zeros(1,N);
    v_m = zeros(1,N);
    iteration = 1;
    while t <= tEnd
       dvdt = zeros(1,N);
       dedt = zeros(1,N);
       if iteration > 1         % leap frog for time advance
           v_m(20:end-20) = v(20:end-20);
           e_m(20:end-20) = e(20:end-20);
           v(20:end-20) = v(20:end-20) + dvdt(20:end-20) * dt / 2;
           e(20:end-20) = e(20:end-20) + dedt(20:end-20) * dt / 2;
       end
       [rho, v, e, p, x, dvdt, dedt] = single_step(m, N, rho, v, e, p, x, h, dvdt, dedt);
       if iteration == 1
           v(20:end-20) = v(20:end-20) + dvdt(20:end-20) * dt / 2;
           e(20:end-20) = e(20:end-20) + dedt(20:end-20) * dt / 2;
           x(20:end-20) = x(20:end-20) + v(20:end-20) * dt;
       else
           v(20:end-20) = v_m(20:end-20) + dvdt(20:end-20) * dt;
           e(20:end-20) = e_m(20:end-20) + dedt(20:end-20) * dt;
           x(20:end-20) = x(20:end-20) + v(20:end-20) * dt;
       end
       h = h0 .* rho0 ./ rho;
       t = t + dt;
       iteration = iteration + 1;
    end
    subplot(2,2,1);plot(x,rho,'o','markersize',3);xlim([-0.6 0.6]);ylabel('Density');xlabel('Position');
    subplot(2,2,2);plot(x,v,'o','markersize',3);xlim([-0.6 0.6]);ylabel('Velocity');xlabel('Position');
    subplot(2,2,3);plot(x,p,'o','markersize',3);xlim([-0.6 0.6]);ylabel('Pressure');xlabel('Position')
    subplot(2,2,4);plot(x,e,'o','markersize',3);xlim([-0.6 0.6]);ylabel('Internal energy');xlabel('Position');
end

function [rho, v, e, p, x, h, dvdt, dedt] = single_step(m, N, rho, v, e, p, x, h, dvdt, dedt)
   global gamma
   [PAIR_I, PAIR_J, NIAC, W, DWDX] = NNPS(x,N,h);
   % compute density
   rho = m * 2/3 ./ h;  % Wii
   for k = 1:NIAC
       i = PAIR_I(k);
       j = PAIR_J(k);
       rho(i) = rho(i) + m*W(k);
       rho(j) = rho(j) + m*W(k);    
   end   
   p = (gamma-1).*rho.*e;     % compute pressure
   c = sqrt(gamma.*(gamma-1)*e);   % compute sound speed
   for k = 1:NIAC
       i = PAIR_I(k);
       j = PAIR_J(k);
       cij = (c(i)+c(j)) / 2;
       hij = (h(i)+h(j)) / 2;
       vij = v(i)-v(j);
       xij = x(i)-x(j);
       phiij = (hij * vij * xij) / (xij^2 + (0.1*hij)^2);
       rhoij = (rho(i)+rho(j))/2;
       if vij*xij < 0
            Piij = p(i)/rho(i)^2 + p(j)/rho(j)^2 + (-0.5 * cij * phiij + phiij^2) / rhoij; % artificial viscosity
       else
            Piij = p(i)/rho(i)^2 + p(j)/rho(j)^2;
       end
       dvdt(i) = dvdt(i) + m * Piij * DWDX(k);
       dvdt(j) = dvdt(j) - m * Piij * DWDX(k);
       dedt(i) = dedt(i) + m * Piij * DWDX(k) * vij / 2;
       dedt(j) = dedt(j) + m * Piij * DWDX(k) * vij / 2;
   end
end

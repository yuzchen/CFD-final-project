function [rho0, u0, p0, tEnd, x0] = IC_setup(x, initType)
switch initType
    case{1} 
        rho = [1.0 0.125];
        u = [0 0];
        p = [1 0.1];
        tEnd = 0.14; 
        x0 = 0;
        disp(['initType = ', num2str(initType), ', rho = ', num2str(rho(1)),' ', num2str(rho(2)),...
            ', u = ', num2str(u(1)),' ', num2str(u(2)), ', p = ', num2str(p(1)),' ', num2str(p(2)), ', tEnd = ', num2str(tEnd), ', x0 = ', num2str(x0)]);     
	case{2} 
        rho = [1 0.125];
        u = [0.75 0];
        p = [1 0.1];
        tEnd = 0.2; 
        x0 = -0.2;
        disp(['initType = ', num2str(initType), ', rho = ', num2str(rho(1)),' ', num2str(rho(2)),...
            ', u = ', num2str(u(1)),' ', num2str(u(2)), ', p = ', num2str(p(1)),' ', num2str(p(2)), ', tEnd = ', num2str(tEnd), ', x0 = ', num2str(x0)]);
    case{3} 
        rho = [0.445 0.5];
        u = [0.698 0];
        p = [3.528 0.571];
        tEnd = 0.1; 
        x0 = 0.5;
        disp(['initType = ', num2str(initType), ', rho = ', num2str(rho(1)),' ', num2str(rho(2)),...
            ', u = ', num2str(u(1)),' ', num2str(u(2)), ', p = ', num2str(p(1)),' ', num2str(p(2)), ', tEnd = ', num2str(tEnd), ', x0 = ', num2str(x0)]);  
        case{4} 
        rho = [1 0.25];
        u = [0 0];
        p = [1 0.1795];
        tEnd = 0.2; 
        x0 = 0;
        disp(['initType = ', num2str(initType), ', rho = ', num2str(rho(1)),' ', num2str(rho(2)),...
            ', u = ', num2str(u(1)),' ', num2str(u(2)), ', p = ', num2str(p(1)),' ', num2str(p(2)), ', tEnd = ', num2str(tEnd), ', x0 = ', num2str(x0)]);    
    otherwise 
        error('wrong initType');
end
left = find(x <= x0);
right = find(x > x0);
rho0(left) = rho(1);
rho0(right) = rho(2);
u0(left) = u(1);
u0(right) = u(2);
p0(left) = p(1);
p0(right) = p(2);
end
function [delta0,E0,TM0] = getGenSteadyStateVal(v0,theta0,pg0,qg0,gen_set,xq,xdp,nG)
% Get generator's steady-state values at equilibrium
% Author: Sebastian A. Nugroho
% Date: 1/27/2021

%Solve the nonlinear function
problem.options = optimoptions('fsolve','Display','off');
problem.objective = @(x)genStateFun(x,v0,theta0,pg0,qg0,gen_set,xq,xdp,nG);
problem.x0 = zeros(2*nG,1); 
problem.solver='fsolve'; 
x = fsolve(problem); 

%Extract the results
delta0 = x(1:2:2*nG);
E0 = x(2:2:2*nG);

%Compute the remaining values: TMi 
TM0 = pg0;

function val = genStateFun(x,v0,theta0,pg0,qg0,gen_set,xq,xdp,nG)
% The variable x is constructed as x = [x1, x2, ..., x_nG]
% Where x_i = [delta0_i, E0_i];

    v0red = v0(gen_set);
    theta0red = theta0(gen_set);
    for i = 1:nG
        idx_i = (i-1)*2;
        val(idx_i+1) = -pg0(i) + (1/xdp(i))*x(idx_i+2)*v0red(i)*sin(x(idx_i+1)-theta0red(i)) ...
                        - ((xq(i)-xdp(i))/(2*xdp(i)*xq(i)))*v0red(i)^2*sin(2*(x(idx_i+1)-theta0red(i)));
        val(idx_i+2) = -qg0(i) + (1/xdp(i))*x(idx_i+2)*v0red(i)*cos(x(idx_i+1)-theta0red(i)) ...
                        - ((xq(i)+xdp(i))/(2*xdp(i)*xq(i)))*v0red(i)^2 ...
                        - ((xq(i)-xdp(i))/(2*xdp(i)*xq(i)))*v0red(i)^2*cos(2*(x(idx_i+1)-theta0red(i)));
    end
end

end
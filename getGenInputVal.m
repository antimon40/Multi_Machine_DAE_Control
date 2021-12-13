function [EFd0,Tr0] = getGenInputVal(v0,theta0,delta0,E0,TM0,gen_set,xd,xdp)
% Get generator's input values at equilibrium
% Author: Sebastian A. Nugroho
% Date: 1/28/2021

EFd0 = (xd./xdp).*E0 - ((xd-xdp)./xdp).*v0(gen_set).*cos(delta0-theta0(gen_set));
Tr0 = TM0;

end
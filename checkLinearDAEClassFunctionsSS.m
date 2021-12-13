function [isLinDAEok] = checkLinearDAEClassFunctionsSS(sys,Xss,Uss)
% Check the linear DAE for the steady state values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumption: generator buses not connected to loads (fixed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using the power balance equations
% Author: Sebastian A. Nugroho
% Date: 8/4/2020

%Generator dynamics test
gen_dyn_err = full(sys.A_sp)*Xss + full(sys.B_sp)*Uss;

if max(gen_dyn_err) < 1e-1
    disp('Nonlinear DAE equations satisfied.'); 
    isLinDAEok = 1;
else
    disp('Nonlinear DAE equations NOT satisfied.'); 
    isLinDAEok = 0;
end

end
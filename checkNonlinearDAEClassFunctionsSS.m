function [isNonlinDAEok,res] = checkNonlinearDAEClassFunctionsSS(case_sys,Xss,Igss,Vgss,Vlss,Uss)
% Check the nonlinear DAE for the steady state values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumption: generator buses not connected to loads (fixed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using the power balance equations
% Author: Sebastian A. Nugroho
% Date: 8/4/2020

%Generator dynamics test
gen_dyn_err = case_sys.f_gen_dyn(Xss,Igss,Uss);

%Generator algebraic equations test
gen_alg_err = case_sys.f_gen_alg(Xss,Igss,Vgss);

%Power balance equations test
pwr_blc_err = case_sys.f_pf_P(Xss,Igss,Vgss,Vlss);

if sum(gen_dyn_err) + sum(gen_alg_err) + sum(pwr_blc_err) < 1e-1
    disp('Nonlinear DAE equations satisfied.'); 
    isNonlinDAEok = 1;
else
    disp('Nonlinear DAE equations NOT satisfied.'); 
    isNonlinDAEok = 0;
end

%residue
res = sum(gen_dyn_err) + sum(gen_alg_err) + sum(pwr_blc_err);

end
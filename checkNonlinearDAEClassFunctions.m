function [isNonlinDAEok] = checkNonlinearDAEClassFunctions(case_sys)
% Check the nonlinear DAE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumption: generator buses not connected to loads (fixed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using the power balance equations
% Author: Sebastian A. Nugroho
% Date: 7/23/2020

%Initialize vectors
X = [];
Ig = [];
U = [];
VRG0 = case_sys.vR0(case_sys.gen_set); 
VIG0 = case_sys.vI0(case_sys.gen_set);
Vg = [];
for i = 1:case_sys.nG
    X = [X; case_sys.deltaG0(i); case_sys.w0v(i); case_sys.EQIp0(i); case_sys.EDIp0(i)];
    Ig = [Ig; real(case_sys.IDQG0(i)); imag(case_sys.IDQG0(i))];
    U = [U; case_sys.TM0(i); case_sys.EFD0(i)];
    Vg = [Vg; VRG0(i); VIG0(i)];
end 

VRL0 = case_sys.vR0(case_sys.load_set); 
VIL0 = case_sys.vI0(case_sys.load_set);
Vl = [];
for i = 1:case_sys.nL
    Vl = [Vl; VRL0(i); VIL0(i)];
end 

%Generator dynamics test
gen_dyn_err = case_sys.f_gen_dyn(X,Ig,U);

%Generator algebraic equations test
gen_alg_err = case_sys.f_gen_alg(X,Ig,Vg);

%Power balance equations test
pwr_blc_err = case_sys.f_pf_P(X,Ig,Vg,Vl);

if max(gen_dyn_err) + sum(gen_alg_err) + sum(pwr_blc_err) < 1e-1
    disp('Nonlinear DAE equations satisfied.'); 
    isNonlinDAEok = 1;
else
    disp('Nonlinear DAE equations NOT satisfied.'); 
    isNonlinDAEok = 0;
end

end
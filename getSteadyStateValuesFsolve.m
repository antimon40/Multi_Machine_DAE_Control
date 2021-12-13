function [Z_ss,X_ss,Ig_ss,Vg_ss,Vl_ss,U] = getSteadyStateValuesFsolve(case_sys)
% Get the steady state values using Newton's method using Fsolve
% Author: Sebastian A. Nugroho
% Date: 8/13/2020

%Jacobian function name
name_jacobian_m = sprintf('Jac_Fsolve_%s.m',case_sys.casename);
name_jacobian = sprintf('Jac_Fsolve_%s',case_sys.casename);

%Construct function handle for Jacobian matrix
jac_handle = str2func(name_jacobian);

%Bus voltages
vR0g = case_sys.vR0(case_sys.gen_set);
vI0g = case_sys.vI0(case_sys.gen_set);
vR0l = case_sys.vR0(case_sys.load_set);
vI0l = case_sys.vI0(case_sys.load_set);

%Create Z0: vector of initial conditions
%Construct X0,Ig0,Vg0,U
X0 = zeros(4*case_sys.nG,1);
Ig0 = zeros(2*case_sys.nG,1);
Vg0 = zeros(2*case_sys.nG,1);
U = zeros(2*case_sys.nG,1);
for i = 1:case_sys.nG
    idx_X = (i-1)*4;
    idx_Ig = (i-1)*2;
    idx_Vg = (i-1)*2;
    idx_U = (i-1)*2;
    X0(idx_X+1) = case_sys.deltaG0(i);
    X0(idx_X+2) = case_sys.w0v(i);
    X0(idx_X+3) = case_sys.EQIp0(i);
    X0(idx_X+4) = case_sys.EDIp0(i);
    Ig0(idx_Ig+1) = real(case_sys.IDQG0(i));
    Ig0(idx_Ig+2) = imag(case_sys.IDQG0(i));
    Vg0(idx_Vg+1) = vR0g(i);
    Vg0(idx_Vg+2) = vI0g(i);
    U(idx_U+1) = case_sys.TM0(i);
    U(idx_U+2) = case_sys.EFD0(i);
end
Vl0 = zeros(2*case_sys.nL,1);
for i = 1:case_sys.nL
    idx_Vl = (i-1)*2;
    Vl0(idx_Vl+1) = vR0l(i);
    Vl0(idx_Vl+2) = vI0l(i);
end
Z0 = [X0; Ig0; Vg0; Vl0];

%Fsolve options
options = optimoptions('fsolve','Algorithm','trust-region','SpecifyObjectiveGradient',true,...
    'UseParallel',false,'FunctionTolerance',1.0000e-12,'MaxIter',10000,'StepTolerance', 1.0000e-12,...
    'MaxFunctionEvaluations',1e4);

%Solve using Fsolve
[Z_ss,fval] = fsolve(@(z)fun_all(z), [Z0], options);

%Show fval
fprintf(1, 'Show fval norm: %f',norm(fval,2));
fprintf('\n');

%Results
X_ss = Z_ss(1:4*case_sys.nG);
Ig_ss = Z_ss(4*case_sys.nG+1:(4+2)*case_sys.nG); 
Vg_ss = Z_ss((4+2)*case_sys.nG+1:(4+2+2)*case_sys.nG); 
Vl_ss = Z_ss((4+2+2)*case_sys.nG+1:end); 

function [fz,jacobianz] = fun_all(z)
% function [fz] = fun_all(z)
    Xs = z(1:4*case_sys.nG);
    Igs = z(4*case_sys.nG+1:(4+2)*case_sys.nG); 
    Vgs = z((4+2)*case_sys.nG+1:(4+2+2)*case_sys.nG); 
    Vls = z((4+2+2)*case_sys.nG+1:end); 
    fz = case_sys.F_all(Xs,Igs,Vgs,Vls,U);
    jacobianz = feval(jac_handle,Xs,Igs,Vgs,Vls); 
    jacobianz = jacobianz';
end

end
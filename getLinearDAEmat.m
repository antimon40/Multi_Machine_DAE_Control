function [Ad_d,Ad_a,Aa_d,Aa_a,Bd,Bq] = getLinearDAEmat(case_sys)
% Get the matrices for the linearized DAE around the equilibrium
% Author: Sebastian A. Nugroho
% Date: 2/18/2021

%Save the Jacobian matrix into a file
name_jacobian1 = sprintf('Jac_Ad_d_%s',case_sys.casename);
name_jacobian2 = sprintf('Jac_Ad_a_%s',case_sys.casename);
name_jacobian3 = sprintf('Jac_Bd_%s',case_sys.casename);
name_jacobian4 = sprintf('Jac_Aa_d_%s',case_sys.casename);
name_jacobian5 = sprintf('Jac_Aa_a_%s',case_sys.casename);
name_jacobian6 = sprintf('Jac_Bq_%s',case_sys.casename);

%Construct function handle for Jacobian matrix
jac_handle1 = str2func(name_jacobian1);
jac_handle2 = str2func(name_jacobian2);
jac_handle3 = str2func(name_jacobian3);
jac_handle4 = str2func(name_jacobian4);
jac_handle5 = str2func(name_jacobian5);
jac_handle6 = str2func(name_jacobian6);

%Auxiliary variable
w0 = case_sys.w0*ones(case_sys.nG,1);

%Compute the matrices 
Ad_d = feval(jac_handle1,case_sys.delta0,w0,case_sys.Ep0,case_sys.TM0,...
    case_sys.pg0,case_sys.qg0,case_sys.v0,case_sys.theta0);
Ad_a = feval(jac_handle2,case_sys.delta0,w0,case_sys.Ep0,case_sys.TM0,...
    case_sys.pg0,case_sys.qg0,case_sys.v0,case_sys.theta0);
Aa_d = feval(jac_handle4,case_sys.delta0,w0,case_sys.Ep0,case_sys.TM0,...
    case_sys.pg0,case_sys.qg0,case_sys.v0,case_sys.theta0);
Aa_a = feval(jac_handle5,case_sys.delta0,w0,case_sys.Ep0,case_sys.TM0,...
    case_sys.pg0,case_sys.qg0,case_sys.v0,case_sys.theta0);
Bd = feval(jac_handle3,case_sys.EFd0,case_sys.Tr0);
Bq = feval(jac_handle6);

end
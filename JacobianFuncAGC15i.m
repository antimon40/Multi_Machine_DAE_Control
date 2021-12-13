function [JacMat,JacMatprime] = JacobianFuncAGC15i(t,state,stateprime)
% Get the Jacobian matrix
% Author: Sebastian A. Nugroho
% Date: 2/21/2021

%Load power network's object
filename = 'PNobj.mat';
load(filename);

%Jacobian function name
name_jacobian = sprintf('Jac_%s',case_sys.casename);
name_jacobian0 = sprintf('Jac_u_%s',case_sys.casename);
name_jacobian1 = sprintf('Jac_Ad_d_%s',case_sys.casename);
name_jacobian2 = sprintf('Jac_Ad_a_%s',case_sys.casename);
name_jacobian4 = sprintf('Jac_Aa_d_%s',case_sys.casename);
name_jacobian5 = sprintf('Jac_Aa_a_%s',case_sys.casename);

%Construct function handle for Jacobian matrix
jac_handle = str2func(name_jacobian);
jac_handle0 = str2func(name_jacobian0);
jac_handle1 = str2func(name_jacobian1);
jac_handle2 = str2func(name_jacobian2);
jac_handle4 = str2func(name_jacobian4);
jac_handle5 = str2func(name_jacobian5);

%Reshape vectors
dels = state(1:case_sys.nG,1);
oms = state(case_sys.nG+1:2*case_sys.nG,1);
Eps = state(2*case_sys.nG+1:3*case_sys.nG,1);
Tms = state(3*case_sys.nG+1:4*case_sys.nG,1);
pgs = state(4*case_sys.nG+1:5*case_sys.nG,1);
qgs = state(5*case_sys.nG+1:6*case_sys.nG,1);
vs = state(6*case_sys.nG+1:6*case_sys.nG+case_sys.nN,1);
thes = state(6*case_sys.nG+1+case_sys.nN:end,1);

%Compute the matrices 
Ad_d = full(feval(jac_handle1,dels,oms,Eps,Tms,pgs,qgs,vs,thes));
Ad_a = full(feval(jac_handle2,dels,oms,Eps,Tms,pgs,qgs,vs,thes));
Aa_d = full(feval(jac_handle4,dels,oms,Eps,Tms,pgs,qgs,vs,thes));
Aa_a = full(feval(jac_handle5,dels,oms,Eps,Tms,pgs,qgs,vs,thes));

%Additional terms due to AGC
dfydy = -case_sys.Kg;
dfydw = -(case_sys.Kg/case_sys.nG)*(1./case_sys.Rd + case_sys.D)';
dfydpg = case_sys.Kg*ones(1,case_sys.nG);
add_terms_row = [zeros(1,case_sys.nG) dfydw zeros(1,2*case_sys.nG) dfydy ...
            dfydpg zeros(1,1*case_sys.nG) zeros(1,2*case_sys.nN)];
add_terms_col1 = zeros(4*case_sys.nG,1);
add_terms_col2 = zeros(2*case_sys.nG+2*case_sys.nN,1);

%Compute the Jacobian matrix 
Jac_full = [Ad_d add_terms_col1 Ad_a; add_terms_row; Aa_d add_terms_col2 Aa_a];
Jac_sparse = sparse(Jac_full);

temp = feval(jac_handle0,case_sys.EFd0,case_sys.Tr0);
temp(end+1,:) = zeros(1,2*case_sys.nG);

Jac_sparse2 = temp*sparse(case_sys.Kagcfull);

JacMat = Jac_sparse + Jac_sparse2;

%Load the mass matrix
filename = 'MassMatrixAGC.mat';
load(filename);

JacMatprime = -MassAGC;

end
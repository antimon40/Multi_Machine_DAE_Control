function [Jac_sparse] = JacobianAlgFunc(state,delta0plus,E0plus,pd0plus,qd0plus)
% Get the Jacobian matrix for the algebraic equations
% Author: Sebastian A. Nugroho
% Date: 1/31/2021

%Load power network's object
filename = 'PNobj.mat';
load(filename);

%Jacobian function name
name_jacobian = sprintf('Jac_alg_%s',case_sys.casename);

%Construct function handle for Jacobian matrix
jac_handle = str2func(name_jacobian);

%Reshape vectors
pgs = state(1:case_sys.nG,1);
qgs = state(case_sys.nG+1:2*case_sys.nG,1);
vs = state(2*case_sys.nG+1:2*case_sys.nG+case_sys.nN,1);
thes = state(2*case_sys.nG+case_sys.nN+1:end,1);

%Compute the Jacobian matrix 
Jac_sparse = feval(jac_handle,delta0plus,E0plus,pgs,qgs,vs,thes);

end
function [Jac_sparse] = JacobianFunc(t,state)
% Get the Jacobian matrix
% Author: Sebastian A. Nugroho
% Date: 1/30/2021

%Load power network's object
filename = 'PNobj.mat';
load(filename);

%Jacobian function name
name_jacobian = sprintf('Jac_%s',case_sys.casename);

%Construct function handle for Jacobian matrix
jac_handle = str2func(name_jacobian);

%Reshape vectors
dels = state(1:case_sys.nG,1);
oms = state(case_sys.nG+1:2*case_sys.nG,1);
Eps = state(2*case_sys.nG+1:3*case_sys.nG,1);
Tms = state(3*case_sys.nG+1:4*case_sys.nG,1);
pgs = state(4*case_sys.nG+1:5*case_sys.nG,1);
qgs = state(5*case_sys.nG+1:6*case_sys.nG,1);
vs = state(6*case_sys.nG+1:6*case_sys.nG+case_sys.nN,1);
thes = state(6*case_sys.nG+1+case_sys.nN:end,1);

%Compute the Jacobian matrix 
Jac_sparse = feval(jac_handle,dels,oms,Eps,Tms,pgs,qgs,vs,thes);

end
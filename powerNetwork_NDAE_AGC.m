function [val] = powerNetwork_NDAE_AGC(t,state,input,Pl,Ql,case_sys)
% Construct the overall power network's NDAE model with AGC
% Author: Sebastian A. Nugroho
% Date: 2/22/2021

%Reshape vectors - state
dels = state(1:case_sys.nG,1);
oms = state(case_sys.nG+1:2*case_sys.nG,1);
Eps = state(2*case_sys.nG+1:3*case_sys.nG,1);
Tms = state(3*case_sys.nG+1:4*case_sys.nG,1);
ys = state(4*case_sys.nG+1,1);
pgs = state(4*case_sys.nG+2:5*case_sys.nG+1,1);
qgs = state(5*case_sys.nG+2:6*case_sys.nG+1,1);
vs = state(6*case_sys.nG+2:6*case_sys.nG+case_sys.nN+1,1);
thes = state(6*case_sys.nG+2+case_sys.nN:end,1);

%Reshape vectors - input
Efd = input(1:case_sys.nG,1);
Tr = input(case_sys.nG+1:end,1);

%Function #1: generator's differential equations
[out] = case_sys.f_gen_dyn_AGC(dels,oms,Eps,Tms,Efd,Tr,pgs,vs,thes,ys);

%Function #2: generator's algebraic equations
[val1,val2] = case_sys.f_gen_alg(dels,Eps,pgs,qgs,vs,thes);

%Function #3: power flow equations
[valGp,valGq,valLp,valLq] = case_sys.f_pf_loads(pgs,qgs,vs,thes,Pl,Ql);

%Output
val = [out; val1; val2; valGp; valGq; valLp; valLq];

end
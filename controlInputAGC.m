function [u] = controlInputAGC(case_sys,sys,initial_input,state,K,x0d)
% Construct control inputs
% Author: Sebastian A. Nugroho
% Date: 2/22/2021

feedback = K*(state(1:sys.nxd+1) - x0d);
Kpart = K(:,4*case_sys.nG+1);

u = [initial_input(1:case_sys.nG,1) + feedback(1:case_sys.nG,1); ...
     initial_input(case_sys.nG+1:end,1) + Kpart(case_sys.nG+1:end,1)*state(4*case_sys.nG+1,1)];

end
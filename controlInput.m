function [u] = controlInput(case_sys,sys,initial_input,state,K)
% Construct control inputs
% Author: Sebastian A. Nugroho
% Date: 1/30/2021

u = initial_input + K*(state(1:sys.nxd) - case_sys.x0);

end
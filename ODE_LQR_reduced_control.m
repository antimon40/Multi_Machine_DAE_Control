function [Kd,K,elapsedTime] = ODE_LQR_reduced_control(sys,x0)
% Computing a stabilizing controller gain matrix for NDAE
% Based on LQR and linearized ODE power networks
% Bring the states to the zero equilibrium
% Author: Sebastian A. Nugroho
% Date: 2/19/2021

%% Controller Design
%Constants for positive definiteness
eps1 = 1e-9; 
eps2 = 1e-9;

%Cost matrices
Q = 1e-1*eye(sys.nxlin);
R = 1e-1*eye(sys.nulin);

%Solve the LMI
%Optimization variables
P = sdpvar(sys.nxlin,sys.nxlin,'symmetric');
gamma = sdpvar(1);

%Objective functions
obj =  gamma;

F1 = [sys.Alin*P + P*sys.Alin' - sys.Bulin*inv(R)*sys.Bulin' P; P -inv(Q)] + eps1*eye(2*sys.nxlin) <= 0;

F2 = [gamma x0'; x0 P] - eps1*eye(sys.nxlin + 1) >= 0;

F3 = [P >= eps2*eye(sys.nxlin), gamma >= eps2];

%YALMIP settings
ops = sdpsettings('verbose',1,'solver','mosek','showprogress',1);

%Solve the optimization problem
disp('Solving LMIs...');
tt = tic;
sol = optimize([F1, F2, F3], obj, ops);
elapsedTime = toc(tt);
disp('Done!');
fprintf(1, 'Computation time: %f seconds', elapsedTime);
fprintf(1, '\nYALMIP status: %s', sol.info);
fprintf('\n');

%Values
P = value(P);
gamma = value(gamma);

%Replace NaN with zeros
P(isnan(P)) = 0;
gamma(isnan(gamma)) = 0;

%Matrices
Kd = -inv(R)*sys.Bulin'*inv(P);
K = [Kd zeros(sys.nud,sys.nxa)];



end
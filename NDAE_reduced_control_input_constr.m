function [Kd,K] = NDAE_reduced_control_input_constr(sys,xd0)
% Computing a stabilizing controller gain matrix for NDAE
% Controller is ud = Kd*xd
% The control input is constrained
% Bring the states to the zero equilibrium
% Author: Sebastian A. Nugroho
% Date: 2/13/2021

%% Controller Design
%Constants for positive definiteness
eps1 = 1e-12; 

%Tuning parameters
mu_bar = 1e+12; %for Case 9
R = 1e-4*eye(sys.nud);
umax = 1e+5;
% mu_bar = 1e+13; %for Case 14
% R = 1e-4*eye(sys.nud);
% umax = 1e+7;

%Solve the LMI
%Optimization variables
X1 = sdpvar(sys.nxd,sys.nxd,'symmetric');
X2 = sdpvar(sys.nxa,sys.nxd);
Y1 = sdpvar(sys.nxa,sys.nxd);
Y2 = sdpvar(sys.nxa,sys.nxa);
W = sdpvar(sys.nud,sys.nxd);
sigma = sdpvar(1);

%Objective functions
obj =  []; %1e-4*(norm(W,2));

F1 = [sys.Ed'*X1*sys.Ad' + sys.Ad*X1*sys.Ed + sys.Ed'*W'*sys.Bud' + sys.Bud*W*sys.Ed + sigma*sys.Gd*sys.Gd' ... %1st row
    sys.Ed'*X2'*sys.Aa' + Y1'*sys.Aa' sys.Ed'*X1*sys.Hd' sys.Ed'*X2'*sys.Ha' + Y1'*sys.Ha' ; ... %1st row
    sys.Aa*X2*sys.Ed + sys.Aa*Y1 Y2'*sys.Aa' + sys.Aa*Y2 + sigma*sys.Ga*sys.Ga' zeros(sys.nxa,sys.nxd) Y2'*sys.Ha' ; ...
    sys.Hd*X1*sys.Ed zeros(sys.nxd,sys.nxa) -sigma*eye(sys.nxd) zeros(sys.nxd,sys.nxa); ...
    sys.Ha*X2*sys.Ed + sys.Ha*Y1 sys.Ha*Y2 zeros(sys.nxa,sys.nxd) -sigma*eye(sys.nxa)] + eps1*eye(2*sys.nxd + 2*sys.nxa) <= 0;

F2 = [X1 >= eps1*eye(sys.nxd), sigma >= eps1];

F3 = [-mu_bar*eye(sys.nxd) eye(sys.nxd); eye(sys.nxd) -X1*sys.Ed*inv(sys.Ed)'] + eps1*eye(2*sys.nxd) <= 0; 

F4 = [-(umax/(mu_bar*norm(xd0,2)))*X1*sys.Ed*inv(sys.Ed)' inv(sys.Ed)*sys.Ed'*W'; W*sys.Ed*inv(sys.Ed)' -inv(R)] + eps1*eye(sys.nxd + sys.nud) <= 0; 
  
%YALMIP settings
ops = sdpsettings('verbose',1,'solver','mosek','showprogress',1);

%Solve the optimization problem
disp('Solving LMIs...');
tt = tic;
sol = optimize([F1, F2, F3, F4], obj, ops);
elapsedTime = toc(tt);
disp('Done!');
fprintf(1, 'Computation time: %f seconds', elapsedTime);
fprintf(1, '\nYALMIP status: %s', sol.info);
fprintf('\n');

%Values
Y1 = value(Y1);
X1 = value(X1);
W = value(W);
sigma = value(sigma);

%Replace NaN with zeros
Y1(isnan(Y1)) = 0;
X1(isnan(X1)) = 0;
W(isnan(W)) = 0;

%Matrices
Kd = W*inv(X1) ;
K = [Kd zeros(sys.nud,sys.nxa)];

end
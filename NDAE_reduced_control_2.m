function [Kd,K] = NDAE_reduced_control_2(sys)
% Computing a stabilizing controller gain matrix for NDAE, version 2
% Controller is ud = Kd*xd
% Bring the states to the zero equilibrium
% Author: Sebastian A. Nugroho
% Date: 3/6/2021

%% Controller Design
%Constants for positive definiteness
eps1 = 1e-13; 
eps2 = 1e-13; 

%Constants
beta1 = 1e-4;
beta2 = 1e-4;

%Solve the LMI
%Optimization variables
X1 = sdpvar(sys.nxd,sys.nxd,'symmetric');
X2 = sdpvar(sys.nxa,sys.nxd);
Y1 = sdpvar(sys.nxa,sys.nxd);
Y2 = sdpvar(sys.nxa,sys.nxa);
W = sdpvar(sys.nud,sys.nxd);
% beta1 = sdpvar(1);
% beta2 = sdpvar(1);
sigma = sdpvar(1);

%Objective functions
obj =  1*1e-4*(norm(W,2)); % + beta1 + beta2;

% F1 = [sys.Ed*X1*sys.Ad' + sys.Ad*X1*sys.Ed' + sys.Ed*W'*sys.Bud' + sys.Bud*W*sys.Ed' + sigma*sys.Gd*sys.Gd' ... %1st row
%     sys.Ed*X2'*sys.Aa' + Y1'*sys.Aa' sys.Ed*X1*sys.Hd' sys.Ed*X2'*sys.Ha' + Y1'*sys.Ha' ; ... %1st row
%     sys.Aa*X2*sys.Ed' + sys.Aa*Y1 Y2'*sys.Aa' + sys.Aa*Y2 + sigma*sys.Ga*sys.Ga' zeros(sys.nxa,sys.nxd) Y2'*sys.Ha' ; ...
%     sys.Hd*X1*sys.Ed' zeros(sys.nxd,sys.nxa) -beta1*sigma*eye(sys.nxd) zeros(sys.nxd,sys.nxa); ...
%     sys.Ha*X2*sys.Ed' + sys.Ha*Y1 sys.Ha*Y2 zeros(sys.nxa,sys.nxd) -beta2*sigma*eye(sys.nxa)] + eps1*eye(2*sys.nxd + 2*sys.nxa) <= 0;

% F1 = [sys.Ed*X1*sys.Ad' + sys.Ad*X1*sys.Ed' + sys.Ed*W'*sys.Bud' + sys.Bud*W*sys.Ed' ... %1st row
%     sys.Ed*X2'*sys.Aa' + Y1'*sys.Aa' sys.Ed*X1*sys.Hd' sys.Ed*X2'*sys.Ha' + Y1'*sys.Ha' sys.Gd zeros(sys.nxd,sys.nga); ... %1st row
%     sys.Aa*X2*sys.Ed' + sys.Aa*Y1 Y2'*sys.Aa' + sys.Aa*Y2 zeros(sys.nxa,sys.nxd) Y2'*sys.Ha' zeros(sys.nxa,sys.ngd) sys.Ga; ...
%     sys.Hd*X1*sys.Ed' zeros(sys.nxd,sys.nxa) -beta1*eye(sys.nxd) zeros(sys.nxd,sys.nxa) zeros(sys.nxd,sys.ngd+sys.nga); ...
%     sys.Ha*X2*sys.Ed' + sys.Ha*Y1 sys.Ha*Y2 zeros(sys.nxa,sys.nxd) -beta2*eye(sys.nxa) zeros(sys.nxa,sys.ngd+sys.nga); ...
%     sys.Gd' zeros(sys.ngd,sys.nxa+sys.nxd+sys.nxa) -eye(sys.ngd) zeros(sys.ngd,sys.nga); ... 
%     zeros(sys.nga,sys.nxd) sys.Ga' zeros(sys.nga,sys.nxd+sys.nxa+sys.ngd) -eye(sys.nga)] + eps1*eye(2*sys.nxd + 2*sys.nxa + sys.ngd + sys.nga) <= 0;

% F1 = [sys.Ed*X1*sys.Ad' + sys.Ad*X1*sys.Ed' + sys.Ed*W'*sys.Bud' + sys.Bud*W*sys.Ed' ... %1st row
%     sys.Ed*X2'*sys.Aa' + Y1'*sys.Aa' sys.Ed*X1*sys.Hd' sys.Ed*X2'*sys.Ha' + Y1'*sys.Ha' sigma*sys.Gd zeros(sys.nxd,sys.nga); ... %1st row
%     sys.Aa*X2*sys.Ed' + sys.Aa*Y1 Y2'*sys.Aa' + sys.Aa*Y2 zeros(sys.nxa,sys.nxd) Y2'*sys.Ha' zeros(sys.nxa,sys.ngd) sigma*sys.Ga; ...
%     sys.Hd*X1*sys.Ed' zeros(sys.nxd,sys.nxa) -beta1*sigma*eye(sys.nxd) zeros(sys.nxd,sys.nxa) zeros(sys.nxd,sys.ngd+sys.nga); ...
%     sys.Ha*X2*sys.Ed' + sys.Ha*Y1 sys.Ha*Y2 zeros(sys.nxa,sys.nxd) -beta2*sigma*eye(sys.nxa) zeros(sys.nxa,sys.ngd+sys.nga); ...
%     sigma*sys.Gd' zeros(sys.ngd,sys.nxa+sys.nxd+sys.nxa) -sigma*eye(sys.ngd) zeros(sys.ngd,sys.nga); ... 
%     zeros(sys.nga,sys.nxd) sigma*sys.Ga' zeros(sys.nga,sys.nxd+sys.nxa+sys.ngd) -sigma*eye(sys.nga)] + eps1*eye(2*sys.nxd + 2*sys.nxa + sys.ngd + sys.nga) <= 0;

F2 = [X1 >= eps2*eye(sys.nxd), sigma >= eps1]; %beta1 >= eps1, beta2 >= eps1];  
  
%YALMIP settings
ops = sdpsettings('verbose',1,'solver','mosek','showprogress',1);

%Solve the optimization problem
disp('Solving LMIs...');
tt = tic;
sol = optimize([F1, F2], obj, ops);
elapsedTime = toc(tt);
disp('Done!');
fprintf(1, 'Computation time: %f seconds', elapsedTime);
fprintf(1, '\nYALMIP status: %s', sol.info);
fprintf('\n');

%Values
Y1 = value(Y1);
X1 = value(X1);
W = value(W);
beta1 = value(beta1);
beta2 = value(beta2);

%Replace NaN with zeros
Y1(isnan(Y1)) = 0;
X1(isnan(X1)) = 0;
W(isnan(W)) = 0;

%Matrices
Kd = W*inv(X1) ;
K = [Kd zeros(sys.nud,sys.nxa)];

%Compute aux matrices
EdX1 = sys.Ed'*inv(sys.Ed)*inv(X1);

end
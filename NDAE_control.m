function [K] = NDAE_control(sys)
% Computing a stabilizing controller gain matrix for NDAE
% Bring the states to the zero equilibrium
% Author: Sebastian A. Nugroho
% Date: 2/2/2021

%% Construct an orthogonal complement
%Find an orthogonal complement of E
sys.E_comp = null(sys.E');

%Resize
sys.E_comp = sys.E_comp(:,1:sys.nx-sys.re)';

%Check if it is an orthogonal complement of E
I_E = sys.E_comp*sys.E_comp';
if (nnz(I_E - eye(size(I_E,1))) == 0) && (nnz(sys.E_comp*sys.E) == 0)
    disp('An orthogonal complement is found.'); 
else
    disp('An orthogonal complement is NOT found.'); 
    return;
end

%% Controller Design
%Constants for positive definiteness
eps1 = 1e-12; %eps; %1e-16;

%Solve the LMI
%Optimization variables
X = sdpvar(sys.nx,sys.nx,'symmetric');
Y = sdpvar(sys.nx-sys.re,sys.nx);
W = sdpvar(sys.nu,sys.nx);
beta = sdpvar(1);

%Objective functions
obj =  1*(norm(X,2) + norm(Y,2) + norm(W,2));

F1 = [sys.A*X*sys.E + sys.E'*X*sys.A' + sys.A*sys.E_comp'*Y + Y'*sys.E_comp*sys.A' + sys.Bu*W + W'*sys.Bu' + beta*sys.G*sys.G' sys.E'*X*sys.H' + Y'*sys.E_comp*sys.H'; ...
      sys.H*X*sys.E + sys.H*sys.E_comp'*Y                                                                                      -beta*eye(sys.nx)] + eps1*eye(2*sys.nx) <= 0;

F2 = [X >= eps1*eye(sys.nx), beta >= eps1];
  
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
Y = value(Y);
X = value(X);
W = value(W);
beta = value(beta);

%Replace NaN with zeros
Y(isnan(Y)) = 0;
X(isnan(Y)) = 0;
W(isnan(Y)) = 0;

%Matrices
Qinv = inv(X*sys.E + sys.E_comp'*Y);
K = W*Qinv ;

end
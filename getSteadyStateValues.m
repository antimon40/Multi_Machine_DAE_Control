function [Z_k,X_k,Ig_k,Vg_k,Vl_k,U] = getSteadyStateValues(case_sys)
% Get the steady state values using Newton's method based on TI method
% Author: Sebastian A. Nugroho
% Date: 7/24/2020

%Jacobian function name
name_jacobian_m = sprintf('Jac_%s.m',case_sys.casename);
name_jacobian = sprintf('Jac_%s',case_sys.casename);

%Construct function handle for Jacobian matrix
jac_handle = str2func(name_jacobian);

%Bus voltages
vR0g = case_sys.vR0(case_sys.gen_set);
vI0g = case_sys.vI0(case_sys.gen_set);
vR0l = case_sys.vR0(case_sys.load_set);
vI0l = case_sys.vI0(case_sys.load_set);

%Create Z0: vector of initial conditions
%Construct X0,Ig0,Vg0,U
X0 = zeros(4*case_sys.nG,1);
Ig0 = zeros(2*case_sys.nG,1);
Vg0 = zeros(2*case_sys.nG,1);
U = zeros(2*case_sys.nG,1);
for i = 1:case_sys.nG
    idx_X = (i-1)*4;
    idx_Ig = (i-1)*2;
    idx_Vg = (i-1)*2;
    idx_U = (i-1)*2;
    X0(idx_X+1) = case_sys.deltaG0(i);
    X0(idx_X+2) = case_sys.w0v(i);
    X0(idx_X+3) = case_sys.EQIp0(i);
    X0(idx_X+4) = case_sys.EDIp0(i);
    Ig0(idx_Ig+1) = real(case_sys.IDQG0(i));
    Ig0(idx_Ig+2) = imag(case_sys.IDQG0(i));
    Vg0(idx_Vg+1) = vR0g(i);
    Vg0(idx_Vg+2) = vI0g(i);
    U(idx_U+1) = case_sys.TM0(i);
    U(idx_U+2) = case_sys.EFD0(i);
end
Vl0 = zeros(2*case_sys.nL,1);
for i = 1:case_sys.nL
    idx_Vl = (i-1)*2;
    Vl0(idx_Vl+1) = vR0l(i);
    Vl0(idx_Vl+2) = vI0l(i);
end
Z0 = [X0; Ig0; Vg0; Vl0];
    
%Check if the Jacobian function file exists. If it does not, exit
if isfile(name_jacobian_m)
    
    %Set tolerance for termination
    Eps_tol = 10^-6;
    
    %Stopping condition
    stop = false;
    
    %Start the Newton's iteration
    k = 1;
    while ~stop
        %Set prior values for Z
        if k == 1
            Z_p = Z0;
            X_p = X0;
            Ig_p = Ig0;
            Vg_p = Vg0;
            Vl_p = Vl0;
            X_p2 = X0;
            Ig_p2 = Ig0;
        end
        
        %Compute the Jacobian matrix at k-1
        J_p_sparse = feval(jac_handle,X_p,Ig_p,Vg_p,Vl_p);
        
        %Evaluate the negative of F at k-1
        min_F_p_sparse = -case_sys.F_P(X_p,Ig_p,Vg_p,Vl_p,X_p2,Ig_p2,U);
        
        %Solve for Delta_Zk with linsolve
        Delta_Zk = linsolve(full(J_p_sparse),min_F_p_sparse);
        
        %Update Z
        Z_k = Z_p + Delta_Zk;
        
        %Convert to variables
        X_k = Z_k(1:4*case_sys.nG);
        Ig_k = Z_k(4*case_sys.nG+1:(4+2)*case_sys.nG); 
        Vg_k = Z_k((4+2)*case_sys.nG+1:(4+2+2)*case_sys.nG); 
        Vl_k = Z_k((4+2+2)*case_sys.nG+1:end); 
        
        %Check stopping condition
        if norm(Delta_Zk,2) < Eps_tol %Stop
            stop = true;
        else
            X_p2 = X_p;
            Ig_p2 = Ig_p;
            Z_p = Z_k;
            X_p = X_k;
            Ig_p = Ig_k;
            Vg_p = Vg_k;
            Vl_p = Vl_k;
            k = k + 1;
        end
    end
else        
    Z_k = Z0;
    X_k = X0;
    Ig_k = Ig0;
    Vg_k = Vg0;
    Vl_k = Vl0;    
end

end
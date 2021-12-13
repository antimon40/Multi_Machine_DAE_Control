function [isJacobianmade] = getJacobianMatrixFsolve(case_sys)
% Get the Jacobian matrix with constant inputs for Newton's method
% Author: Sebastian A. Nugroho
% Date: 8/13/2020

%Jacobian function name
name_jacobian_m = sprintf('Jac_Fsolve_%s.m',case_sys.casename);

%Check if the Jacobian function file exists. If it does not, build one
if ~isfile(name_jacobian_m)
    %Declare symbolic variables
    Xs = sym('Xs', [4*case_sys.nG 1], 'real');
    Igs = sym('Igs', [2*case_sys.nG 1], 'real');
    Vgs = sym('Vgs', [2*case_sys.nG 1], 'real');
    Vls = sym('Vls', [2*case_sys.nL 1], 'real');

    %Construct U for constant inputs
    U = zeros(2*case_sys.nG,1);
    for i = 1:case_sys.nG
        idx_U = (i-1)*2;
        U(idx_U+1) = case_sys.TM0(i);
        U(idx_U+2) = case_sys.EFD0(i);
    end

    %Redefine function
    fx = case_sys.F_all(Xs,Igs,Vgs,Vls,U);

    %Compute the Jacobian matrix 
    Jmx1 = jacobian(fx,Xs);
    Jmx2 = jacobian(fx,Igs);
    Jmx3 = jacobian(fx,Vgs);
    Jmx4 = jacobian(fx,Vls);
    Jmx = [Jmx1 Jmx2 Jmx3 Jmx4];

    %Save the Jacobian matrix into a file
    name_jacobian = sprintf('Jac_Fsolve_%s',case_sys.casename);
    Jmx_handler = matlabFunction(Jmx,'File',name_jacobian,'Sparse',true,'Vars',{Xs,Igs,Vgs,Vls});
    
    %Output
    disp('Jacobian matrix does not exist.'); 
    disp('A new one has been created.'); 
    isJacobianmade = 1;
else
    Jmx_handler = [];
    name_jacobian = sprintf('Jac_Fsolve_%s',case_sys.casename);
    
    %Output
    disp('Jacobian matrix exists.'); 
    isJacobianmade = 0;
end


end

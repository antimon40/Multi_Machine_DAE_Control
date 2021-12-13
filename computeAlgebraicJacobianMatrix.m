function [isJacobianMade] = computeAlgebraicJacobianMatrix(case_sys)
% Get the Jacobian matrix with constant inputs for the algebraic equations 
% Author: Sebastian A. Nugroho
% Date: 1/31/2021

%Jacobian function name
name_Jacobian_m = sprintf('Jac_alg_%s.m',case_sys.casename);

%Check if the Jacobian function file exists. If it does not, build one
if ~isfile(name_Jacobian_m)
    %Construct the Jacobian matrix
    
    %Declare symbolic variables
    dels = sym('dels', [case_sys.nG 1], 'real'); %delta
%     oms = sym('oms', [case_sys.nG 1], 'real'); %omega
    Eps = sym('Eps', [case_sys.nG 1], 'real'); %Eprime
%     Tms = sym('Tms', [case_sys.nG 1], 'real'); %Tm
    pgs = sym('pgs', [case_sys.nG 1], 'real'); %Pg
    qgs = sym('qgs', [case_sys.nG 1], 'real'); %Qg
    vs = sym('vs', [case_sys.nN 1], 'real'); %v
    thes = sym('thes', [case_sys.nN 1], 'real'); %theta
    
    %Create auxiliary functions
    fgenalg = case_sys.f_gen_alg_compact(dels,Eps,pgs,qgs,vs,thes);
    fpf = case_sys.f_pf_loads_compact(pgs,qgs,vs,thes,case_sys.pl0,case_sys.ql0);
    
    %Start computing the Jacobians
%     dfg2ddel = jacobian(fgenalg,dels);
%     dfg2doms = jacobian(fgenalg,oms);
%     dfg2deps = jacobian(fgenalg,Eps);
%     dfg2dtms = jacobian(fgenalg,Tms);
    dfg2dpgs = jacobian(fgenalg,pgs);
    dfg2dqgs = jacobian(fgenalg,qgs);
    dfg2dvs = jacobian(fgenalg,vs);
    dfg2dthes = jacobian(fgenalg,thes);
%     dpfddel = jacobian(fpf,dels);
%     dpfdoms = jacobian(fpf,oms);
%     dpfdeps = jacobian(fpf,Eps);
%     dpfdtms = jacobian(fpf,Tms);
    dpfdpgs = jacobian(fpf,pgs);
    dpfdqgs = jacobian(fpf,qgs);
    dpfdvs = jacobian(fpf,vs);
    dpfdthes = jacobian(fpf,thes);
    
    %Construct the overall Jacobian matrix of the algebraic equations with respect to the states
    Jmx = [dfg2dpgs dfg2dqgs dfg2dvs dfg2dthes; ...
           dpfdpgs dpfdqgs dpfdvs dpfdthes];
       
    %Save the Jacobian matrix into a file
    name_jacobian = sprintf('Jac_alg_%s',case_sys.casename);
    Jmx_handler = matlabFunction(Jmx,'File',name_jacobian,'Sparse',true,...
        'Vars',{dels,Eps,pgs,qgs,vs,thes});
    
    %Output
    disp('Jacobian matrix does not exist.'); 
    disp('A new one has been created.'); 
    isJacobianMade = 0;
else
    Jmx_handler = [];
    name_jacobian = sprintf('Jac_alg_%s',case_sys.casename);
    
    %Output
    disp('Jacobian matrix exists.'); 
    isJacobianMade = 1;
end

end
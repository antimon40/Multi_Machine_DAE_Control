function [isJacobianMade] = computeJacobianMatrix(case_sys,alwaysCreate)
% Get the Jacobian matrix with constant inputs 
% Author: Sebastian A. Nugroho
% Date: 1/29/2021

%Jacobian function name
name_Jacobian_m = sprintf('Jac_%s.m',case_sys.casename);
name_Jacobian_m2 = sprintf('Jac_u_%s.m',case_sys.casename);
name_Jacobian_m3 = sprintf('Jac_alg_%s.m',case_sys.casename);

%Check if the Jacobian function file exists. If it does not, build one
if (~isfile(name_Jacobian_m) && ~isfile(name_Jacobian_m2) && ~isfile(name_Jacobian_m3)) || alwaysCreate
    %Construct the Jacobian matrix
    
    %Declare symbolic variables
    dels = sym('dels', [case_sys.nG 1], 'real'); %delta
    oms = sym('oms', [case_sys.nG 1], 'real'); %omega
    Eps = sym('Eps', [case_sys.nG 1], 'real'); %Eprime
    Tms = sym('Tms', [case_sys.nG 1], 'real'); %Tm
    pgs = sym('pgs', [case_sys.nG 1], 'real'); %Pg
    qgs = sym('qgs', [case_sys.nG 1], 'real'); %Qg
    vs = sym('vs', [case_sys.nN 1], 'real'); %v
    thes = sym('thes', [case_sys.nN 1], 'real'); %theta
    Efds = sym('Efds', [case_sys.nG 1], 'real'); %Efd
    Trs = sym('Trs', [case_sys.nG 1], 'real'); %Tr
    pds = sym('prs', [case_sys.nN 1], 'real'); %Pd
    qds = sym('qrs', [case_sys.nN 1], 'real'); %Qd
%     ys = sym('ys', [1 1], 'real'); %y
    
    %Create auxiliary functions
    fgendyn = case_sys.f_gen_dyn_compact(dels,oms,Eps,Tms,Efds,Trs,pgs,vs,thes);
%     fgendynAGC = case_sys.f_gen_dyn_AGC(dels,oms,Eps,Tms,Efds,Trs,pgs,vs,thes,ys);
    fgenalg = case_sys.f_gen_alg_compact(dels,Eps,pgs,qgs,vs,thes);
    fpf = case_sys.f_pf_loads_compact(pgs,qgs,vs,thes,pds,qds);
    
    %Start computing the Jacobians with respect to the states
    dfg1ddel = jacobian(fgendyn,dels);
    dfg1doms = jacobian(fgendyn,oms);
    dfg1deps = jacobian(fgendyn,Eps);
    dfg1dtms = jacobian(fgendyn,Tms);
    dfg1dpgs = jacobian(fgendyn,pgs);
    dfg1dqgs = jacobian(fgendyn,qgs);
    dfg1dvs = jacobian(fgendyn,vs);
    dfg1dthes = jacobian(fgendyn,thes);
%     dfg1dys = jacobian(fgendynAGC(end,1),ys);
    dfg2ddel = jacobian(fgenalg,dels);
    dfg2doms = jacobian(fgenalg,oms);
    dfg2deps = jacobian(fgenalg,Eps);
    dfg2dtms = jacobian(fgenalg,Tms);
    dfg2dpgs = jacobian(fgenalg,pgs);
    dfg2dqgs = jacobian(fgenalg,qgs);
    dfg2dvs = jacobian(fgenalg,vs);
    dfg2dthes = jacobian(fgenalg,thes);
    dfg2dpds = jacobian(fgenalg,pds);
    dfg2dqds = jacobian(fgenalg,qds);
    dpfddel = jacobian(fpf,dels);
    dpfdoms = jacobian(fpf,oms);
    dpfdeps = jacobian(fpf,Eps);
    dpfdtms = jacobian(fpf,Tms);
    dpfdpgs = jacobian(fpf,pgs);
    dpfdqgs = jacobian(fpf,qgs);
    dpfdpds = jacobian(fpf,pds);
    dpfdqds = jacobian(fpf,qds);
    dpfdvs = jacobian(fpf,vs);
    dpfdthes = jacobian(fpf,thes);
    
    %Jacobians for dynamics with AGC
    
    %Start computing the Jacobians with respect to the inputs
    dfg1defd = jacobian(fgendyn,Efds);
    dfg1dtr = jacobian(fgendyn,Trs);

    %Construct the overall Jacobian matrix of the NDAE with respect to the states
    Jmx = [dfg1ddel dfg1doms dfg1deps dfg1dtms dfg1dpgs dfg1dqgs dfg1dvs dfg1dthes; ...
           dfg2ddel dfg2doms dfg2deps dfg2dtms dfg2dpgs dfg2dqgs dfg2dvs dfg2dthes; ...
           dpfddel dpfdoms dpfdeps dpfdtms dpfdpgs dpfdqgs dpfdvs dpfdthes];
       
    %Construct the Jacobian matrix of generator's dynamics with respect to the inputs
    Jmu = [dfg1defd dfg1dtr; zeros(2*case_sys.nG,2*case_sys.nG); zeros(2*case_sys.nN,2*case_sys.nG)];
       
    %Construct the overall Jacobian matrix of the algebraic equations with respect to the states
    Jmx2 = [dfg2dpgs dfg2dqgs dfg2dvs dfg2dthes; ...
           dpfdpgs dpfdqgs dpfdvs dpfdthes];
       
    %Compute Jacobians for the linearized DAE model
    Ad_d = [dfg1ddel dfg1doms dfg1deps dfg1dtms];
    Ad_a = [dfg1dpgs dfg1dqgs dfg1dvs dfg1dthes];
    Bd = [dfg1defd dfg1dtr];
    Aa_d = [dfg2ddel dfg2doms dfg2deps dfg2dtms; dpfddel dpfdoms dpfdeps dpfdtms];
    Aa_a = [dfg2dpgs dfg2dqgs dfg2dvs dfg2dthes; dpfdpgs dpfdqgs dpfdvs dpfdthes];
    Bq = [dfg2dpds dfg2dqds; dpfdpds dpfdqds];
       
    %Save the Jacobian matrix into a file
    name_jacobian = sprintf('Jac_%s',case_sys.casename);
    matlabFunction(Jmx,'File',name_jacobian,'Sparse',true,...
        'Vars',{dels,oms,Eps,Tms,pgs,qgs,vs,thes});
    
    %Save the Jacobian matrix into a file
    name_jacobian = sprintf('Jac_u_%s',case_sys.casename);
    matlabFunction(Jmu,'File',name_jacobian,'Sparse',true,...
        'Vars',{Efds,Trs});
    
    %Save the Jacobian matrix into a file
    name_jacobian = sprintf('Jac_alg_%s',case_sys.casename);
    matlabFunction(Jmx2,'File',name_jacobian,'Sparse',true,...
        'Vars',{dels,Eps,pgs,qgs,vs,thes});
    
    %Save the Jacobian matrix into a file
    name_jacobian = sprintf('Jac_Ad_d_%s',case_sys.casename);
    matlabFunction(Ad_d,'File',name_jacobian,'Sparse',true,...
        'Vars',{dels,oms,Eps,Tms,pgs,qgs,vs,thes});
    
    %Save the Jacobian matrix into a file
    name_jacobian = sprintf('Jac_Ad_a_%s',case_sys.casename);
    matlabFunction(Ad_a,'File',name_jacobian,'Sparse',true,...
        'Vars',{dels,oms,Eps,Tms,pgs,qgs,vs,thes});
    
    %Save the Jacobian matrix into a file
    name_jacobian = sprintf('Jac_Bd_%s',case_sys.casename);
    matlabFunction(Bd,'File',name_jacobian,'Sparse',true,...
        'Vars',{Efds,Trs});
    
    %Save the Jacobian matrix into a file
    name_jacobian = sprintf('Jac_Aa_d_%s',case_sys.casename);
    matlabFunction(Aa_d,'File',name_jacobian,'Sparse',true,...
        'Vars',{dels,oms,Eps,Tms,pgs,qgs,vs,thes});
    
    %Save the Jacobian matrix into a file
    name_jacobian = sprintf('Jac_Aa_a_%s',case_sys.casename);
    matlabFunction(Aa_a,'File',name_jacobian,'Sparse',true,...
        'Vars',{dels,oms,Eps,Tms,pgs,qgs,vs,thes});
    
    %Save the Jacobian matrix into a file
    name_jacobian = sprintf('Jac_Bq_%s',case_sys.casename);
    matlabFunction(Bq,'File',name_jacobian,'Sparse',true,...
        'Vars',{});
    
    %Output
    disp('Jacobian matrix does not exist.'); 
    disp('A new one has been created.'); 
    isJacobianMade = 0;
else
    Jmx_handler = [];
    name_jacobian = sprintf('Jac_%s',case_sys.casename);
    
    %Output
    disp('Jacobian matrix exists.'); 
    isJacobianMade = 1;
end

end
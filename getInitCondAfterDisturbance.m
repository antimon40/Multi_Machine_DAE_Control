function [pg0plus,qg0plus,v0plus,theta0plus] = getInitCondAfterDisturbance(case_sys,delta0plus,E0plus,pd0plus,qd0plus)
% Get initial conditions after disturbance occurs
% Author: Sebastian A. Nugroho
% Date: 1/31/2021

options = optimoptions('fsolve','Display','Iter','Algorithm','levenberg-marquardt','InitDamping',0.5, 'ScaleProblem','jacobian',...
    'SpecifyObjectiveGradient',true,'MaxIterations',150,'MaxFunctionEvaluations',300,'OptimalityTolerance',1e-6);

[res,fval,flag] = fsolve(@(x)genStateFun(x,delta0plus,E0plus,pd0plus,qd0plus),...
                [case_sys.pg0; case_sys.qg0; case_sys.v0; case_sys.theta0],options);
            
pg0plus = res(1:case_sys.nG,1);
qg0plus = res(case_sys.nG+1:2*case_sys.nG,1);
v0plus = res(2*case_sys.nG+1:2*case_sys.nG+case_sys.nN,1);
theta0plus = res(2*case_sys.nG+case_sys.nN+1:end,1);
            
if ((flag == 1) || (flag == 2) || (flag == 3) || (flag == 4)) && abs(sum(fval)) < 1e-6
    disp('Initial conditions after disturbance found.'); 
else
    disp('Initial conditions after disturbance NOT found.'); 
end 

function [val,JacVal] = genStateFun(x,delta0plus,E0plus,pd0plus,qd0plus)
    %The variable x is constructed as x = [pg, qg, v, theta]
    
    Pg = x(1:case_sys.nG,1);
    Qg = x(case_sys.nG+1:2*case_sys.nG,1);
    v = x(2*case_sys.nG+1:2*case_sys.nG+case_sys.nN,1);
    theta = x(2*case_sys.nG+case_sys.nN+1:end,1);
    
    [out1] = case_sys.f_gen_alg_compact(delta0plus,E0plus,Pg,Qg,v,theta);
    
    [out2] = case_sys.f_pf_loads_compact(Pg,Qg,v,theta,pd0plus,qd0plus);
    
    val = [out1; out2];
    
    [JacVal] = JacobianAlgFunc(x,delta0plus,E0plus,pd0plus,qd0plus);
    
end

end
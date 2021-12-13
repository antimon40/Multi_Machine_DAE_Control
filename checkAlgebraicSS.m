function [isAlgSSok] = checkAlgebraicSS(case_sys,pd0,qd0)
% Check the steady-state conditions on all algebraic equations
% Author: Sebastian A. Nugroho
% Date: 1/29/2020

%Evaluate function
[val1,val2] = case_sys.f_gen_alg(case_sys.delta0,case_sys.Ep0,case_sys.pg0,case_sys.qg0,case_sys.v0,case_sys.theta0);

[valGp,valGq,valLp,valLq] = case_sys.f_pf(case_sys.pg0,case_sys.qg0,case_sys.v0,case_sys.theta0);
                                   
if abs(sum([sum([val1 val2]) sum([valGp+pd0(case_sys.gen_set) valGq+qd0(case_sys.gen_set)]) ...
        sum([valLp+pd0(case_sys.load_set) valLq+qd0(case_sys.load_set)])])) < 1e-6
    disp('Algebraic steady-state equations satisfied.'); 
    isAlgSSok = 0;
else
    disp('Algebraic steady-state equations NOT satisfied.'); 
    isAlgSSok = 1;
end

end
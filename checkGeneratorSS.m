function [isGenSSok] = checkGeneratorSS(case_sys)
% Check the steady-state conditions on generators
% Author: Sebastian A. Nugroho
% Date: 1/29/2020

%Evaluate function
[delta_dot,omega_dot,Ep_dot,Tm_dot] = case_sys.f_gen_dyn(case_sys.delta0,case_sys.w0*ones(case_sys.nG,1),...
                                        case_sys.Ep0,case_sys.TM0,case_sys.EFd0,case_sys.Tr0,case_sys.pg0,...
                                        case_sys.v0,case_sys.theta0);
                                    
if abs(sum([delta_dot + omega_dot + Ep_dot + Tm_dot])) < 1e-6
    disp('Generator steady-state equations satisfied.'); 
    isGenSSok = 0;
else
    disp('Generator steady-state equations NOT satisfied.'); 
    isGenSSok = 1;
end

end
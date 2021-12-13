function [isPFok] = checkPowerFlowSolution(v0,theta0,pg0,qg0,pl0,ql0,gen_set,load_set,Ybus,nN,nG,nL)
% Check the solution of power flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumption: generator buses not connected to loads (fixed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using the power balance equations
% Author: Sebastian A. Nugroho
% Date: 7/22/2020

%Initialize value
Gdiff = 0;
Ldiff = 0;

%Compute
for i = 1:nN
    if ismember(i,gen_set) %Compute generator buses
        temp = 0;
        for j = 1:nN
            temp = temp + v0(i)*v0(j)*abs(full(Ybus(i,j)))*exp(1i*(theta0(i) - theta0(j) - angle(full(Ybus(i,j)))));
        end
        idx = find(gen_set == i);
        Gdiff = Gdiff + (pg0(idx)+1i*qg0(idx)-pl0(i)-1i*ql0(i)-temp);
    elseif ismember(i,load_set) %Compute load buses
        temp = 0;
        for j = 1:nN
            temp = temp + v0(i)*v0(j)*abs(full(Ybus(i,j)))*exp(1i*(theta0(i) - theta0(j) - angle(full(Ybus(i,j)))));
        end
        Ldiff = Ldiff + (-pl0(i)-1i*ql0(i)-temp);
    end
end

if abs(real(Ldiff+Gdiff) + imag(Ldiff+Gdiff)) < 1e-6
    disp('Power flow equations satisfied.'); 
    isPFok = 0;
else
    disp('Power flow equations NOT satisfied.'); 
    isPFok = 1;
end

end
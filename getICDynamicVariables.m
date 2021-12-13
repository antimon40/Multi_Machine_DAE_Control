function [iG,deltaG,IDQG,VDQG,EDIp,EQIp,EFD,TM,W0] = getICDynamicVariables(v0,theta0,pg0,qg0,pl0,ql0,gen_set,Rs,xd,xq,xdp,xqp,w0,nG)
% Get power system's dynamic variables initial condition values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumption: generator buses not connected to loads (fixed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% According to the steps described in Sauer's book "Power System Dynamics
% and Stability" Section 7.6.3
% Author: Sebastian A. Nugroho
% Date: 7/22/2020

%Initialization
v0gen = v0(gen_set);
theta0gen = theta0(gen_set);
pg0gen = pg0;
qg0gen = qg0;
pl0gen = pl0(gen_set);
ql0gen = ql0(gen_set);

%Step 1: get generator current
iG = (pg0gen + pl0gen -1i.*(qg0gen + ql0gen))./conj((v0gen.*exp(1i*theta0gen)));

%Step 2: get generator current angle delta
deltaG = angle(v0gen.*exp(1i*theta0gen) + (Rs + 1i*xq).*iG);

%Step 3: get generator current and voltage
IDQG = abs(iG).*exp(1i*(angle(iG) - deltaG + (pi/2)*ones(nG,1)));
VDQG = v0gen.*exp(1i*(theta0gen - deltaG + (pi/2)*ones(nG,1)));

%Step 4: compute Ed_prime
EDIp1 = real(VDQG) + Rs.*real(IDQG) - xqp.*imag(IDQG);
EDIp2 = (xq-xqp).*imag(IDQG);

%Check for discrepancy
if (sum(EDIp1-EDIp2)) < 1e-6
    EDIp = EDIp1;
else
    EDIp = (EDIp1 + EDIp2)./2;
end

%Step 5: compute Eq_prime
EQIp = imag(VDQG) + Rs.*imag(IDQG) + xdp.*real(IDQG);

%Step 6: compute Efd
EFD = EQIp + (xd - xdp).*real(IDQG);

%Step 7: compute Tm
TM = EDIp.*real(IDQG) + EQIp.*imag(IDQG) + (xqp - xdp).*real(IDQG).*imag(IDQG);

%Synchronous frequency
W0 = w0*ones(nG,1);

end
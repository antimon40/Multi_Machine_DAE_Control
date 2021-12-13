% Main program - Case 9
% Author: Sebastian A. Nugroho
% Date: 1/25/2020

clear
close all

%Add necessary paths
addpath(genpath('C:\Users\antim\Dropbox\DAE Power Systems\Code\MATPOWER\matpower-master\matpower-master'))
addpath(genpath('C:\Users\antim\Dropbox\DAE Power Systems\Code\PST'))
addpath(genpath('C:\Users\antim\Documents\YALMIP-master\YALMIP-master'))
addpath(genpath('C:\Users\lzg572\Dropbox\MATPOWER\matpower-master\matpower-master'))
addpath(genpath('C:\Users\lzg572\Dropbox\PST'))
addpath(genpath('C:\Users\antim\Dropbox\MATPOWER\matpower-master\matpower-master'))
addpath(genpath('C:\Users\antim\Dropbox\PST'))

%Filename
filename2 = 'Case9matrixDAE.mat';
filename2raw = 'Case9matrixDAE';
casenameraw = 'Case9';

%% Simulation settings

%Need to simulate power systems at steady-state?
isSimulSS = 0;

%Need to simulate power systems after disturbance?
isSimulDist15s = 0;

%Need to simulate power systems after disturbance?
isSimulDist15i = 1;

%Need to stabilize power systems after disturbance?
isRunDist15s = 0;

%Need to stabilize power systems after disturbance?
isRunDist15i = 1;

%Need to stabilize power systems after disturbance?
isRunDistLQR15i = 1;

%Simulation time
tspan1 = 0:0.01:100;

%% Specify renewable and load characteristics
RenMean = 0.2; %20 percent
% RenVariance = 0.01;
RenDisturbance = -0.04; %Low 
LoadDisturbance = 0.04; %Low
distChar = 'low';
% RenDisturbance = -0.08; %Moderate
% LoadDisturbance = 0.08; %Moderate
% distChar = 'moderate';
% RenDisturbance = -0.12; %High
% LoadDisturbance = 0.12; %High
% distChar = 'high';
% RenDisturbance = -0.4; %Ultra
% LoadDisturbance = 0.4; %Ultra
% distChar = 'ultra';

%% Set power system test case
%for MATPOWER
CaseFileNameMATPOWER = 'case9'; %3 machine 9 bus system
%for PST
CaseFileNamePST = 'data3m9b_seb.m'; %3 machine 9 bus system

PauseTime = 1; 

%% Get steady state data from MATPOWER
disp('Get data from MATPOWER.')
%loadcase
network_desc_matpower = loadcase(CaseFileNameMATPOWER);

%Get steady state network description from MATPOWER
[nN,nG,nL,Ybus,Gbus,Bbus,node_set,gen_set,load_set,...
    Ysh,Cf,Ct,Cg,Yff,Yft,Ytf,Ytt,Sbase,Yf,Yt] = getNetworkDescriptionMATPOWER(network_desc_matpower);

%Define index mapping
GenSetMap = [gen_set (1:nG)'];
LoadSetMap = [load_set (1:nL)'];

%% Get generator parameters from PST
disp('Get data from PST.')
%Get generator parameters from PST
[w0,M,D,xq,xd,xqp,xdp,Tq0p,Td0p,Rs,Sbase2] = getNetworkDescriptionPST(CaseFileNamePST);

%Chest valve time constants
% Tch = 0.23*ones(nG,1); %works
% Tch = 0.21*ones(nG,1);
Tch = 0.2*ones(nG,1);

%Regulation constants
% Rd = 0.07*ones(nG,1).*(2*pi); %works
% Rd = 0.07*ones(nG,1).*(2*pi);
Rd = 0.02*ones(nG,1).*(2*pi);

%% Specify the disturbances on renewables and loads
%Assigning 20 % of load at each node to wind generation
%Reducing the overall generation by 20 %
disp('Assigning perturbations to renewables.'); 
ren_set = find(or(network_desc_matpower.bus(:,3)>0, network_desc_matpower.bus(:,4)>0));
pl0 = network_desc_matpower.bus(:,3)./Sbase; % third column of bus matrix is real power demand
ql0 = network_desc_matpower.bus(:,4)./Sbase; % fourth column of bus matrix is reactive power demand
pr0 = zeros(size(pl0));
qr0 = zeros(size(ql0));
pr0(ren_set) = RenMean.*pl0(ren_set);
pd0 = pl0 - pr0;
qd0 = ql0 - qr0;
network_desc_matpower.bus(:,3) = (pd0).*Sbase;
network_desc_matpower.bus(:,4) = (qd0).*Sbase;

%Number of renewables
nR = length(ren_set);

%Nonzero loads
load_set_nonz = ren_set;
nLN = length(load_set_nonz);

%Matrices for mapping renewables and loads
Cr = sparse(nN,nR); 
Cr(sub2ind([nN nR], ren_set, (1:nR).')) = 1;
Cl = sparse(nN,nLN); 
Cl(sub2ind([nN nLN], load_set_nonz, (1:nLN).')) = 1;

%% Get initial values for the state and input variables

%MATPOWER settings for solving power flow
MatpowerOptions = mpoption('model', 'AC', 'pf.alg', 'NR', 'verbose', 0, 'out.all',1,'pf.enforce_q_lims',1,'pf.tol',1e-9); 

disp('Solve the power flow to get initial condition and steady state values.');

%Run the power flow solver via MATPOWER
pf_result = runpf(network_desc_matpower,MatpowerOptions);

%Get initial conditions data from MATPOWER
v0 = pf_result.bus(:,8); % eighth column of bus matrix is voltage magnitude solution
theta0 = degrees2radians(pf_result.bus(:,9)); % nineth column of bus matrix is voltage phase solution
pg0 = pf_result.gen(:,2)./Sbase; % second column of gen matrix is real power set points (or solution for slack bus)
qg0 = pf_result.gen(:,3)./Sbase; % third column of gen matrix is reactive power solutions

%Get difference power
% pd0 = network_desc_matpower.bus(:,3)./Sbase; % third column of bus matrix is real power demand
% qd0 = network_desc_matpower.bus(:,4)./Sbase; % fourth column of bus matrix is reactive power demand

%Sanity check #1
%Check power flow solutions from MATPOWER
[isPFok] = checkPowerFlowSolution(v0,theta0,pg0,qg0,pd0,qd0,gen_set,load_set,Ybus,nN,nG,nL);
pause(PauseTime);

%% Get generator's steady state values
disp('Get generator steady state values.')
[delta0,E0,TM0] = getGenSteadyStateVal(v0,theta0,pg0,qg0,gen_set,xq,xdp,nG);

%% Get generator's inputs at steady state
disp('Get generator inputs.')
[EFd0,Tr0] = getGenInputVal(v0,theta0,delta0,E0,TM0,gen_set,xd,xdp);

%% Assign initial conditions, states, and inputs (at steady-state)
disp('Assign initial conditions.')

%Select order
order_type = 0; 

if order_type == 1  %State-space structure 
    x0 = zeros(4*nG,1);
    a0 = zeros(2*nG,1);
    u0 = zeros(2*nG,1);
    for i = 1:nG
        idx_i = (i-1)*2;
        idx_iG = (i-1)*4;
        x0(idx_iG+1) = delta0(i);
        x0(idx_iG+2) = w0;
        x0(idx_iG+3) = E0(i);
        x0(idx_iG+4) = TM0(i);
        u0(idx_i+1) = EFd0(i);
        u0(idx_i+2) = Tr0(i);
        a0(idx_i+1) = pg0(i);
        a0(idx_i+2) = qg0(i);
    end
    V0 = zeros(2*nN,1);
    V0(1:2:2*nN) = v0;
    V0(2:2:2*nN) = theta0;   
else %Simpler structure
    x0 = [delta0; w0*ones(nG,1); E0; TM0];
    a0 = [pg0; qg0];
    u0 = [EFd0; Tr0];
    V0 = [v0; theta0];    
end

%% Create objects for all generators
disp('Create objects for all generators.')

Genset = {};
if length(M) == nG
    for i = 1:nG
        %Convert generator parameters with the right base MVA (please check if
        %this is correct)
        M_t = M(i)*Sbase2(i)/Sbase;
        D_t = D(i)*Sbase/Sbase2(i);
        xq_t = xq(i)*Sbase/Sbase2(i);
        xd_t = xd(i)*Sbase/Sbase2(i);
        xqp_t = xqp(i)*Sbase/Sbase2(i);
        xdp_t = xdp(i)*Sbase/Sbase2(i);
        Rs_t = Rs(i)*Sbase/Sbase2(i);
        Genset{i} = psgen_class(w0,M_t,D_t,xq_t,xd_t,xqp_t,xdp_t,Tq0p(i),Td0p(i),...
            Tch(i),Rd(i),Rs_t,GenSetMap(i,1),GenSetMap(i,2),Sbase2(i));
    end
end

%% Create a nonlinear DAE system object
disp('Create a nonlinear DAE system object.')

case_sys = nonlin_ps_DAE_class(casenameraw,nN,nG,nL,[],Ybus,Gbus,Bbus,Yf,Yt,Cg,node_set,M,D,xq,xd,xqp,xdp,Tq0p,Td0p,Rs,Tch,Rd,...
                gen_set,load_set,Sbase,w0,Genset,v0,theta0,pg0,qg0,pl0,ql0,pd0,qd0,pr0,qr0,GenSetMap,LoadSetMap,isPFok,delta0,E0,TM0,EFd0,Tr0,x0,u0,a0,V0,...
                Cr,Cl,ren_set,load_set_nonz,nR,nLN);
            
%Save object
filename = 'PNobj.mat';
save(filename,'case_sys');
            
%Sanity check #2
%Check generator's steady-state solutions
[isGenSSok] = checkGeneratorSS(case_sys);
pause(PauseTime);

%Sanity check #3
%Check algebraic constraints steady-state solutions
[isAlgSSok] = checkAlgebraicSS(case_sys,pd0,qd0);
pause(PauseTime);

%% Setting Up AGC
%AGC constants
Kg = 10^3; %Integrator constant
Kpart = pg0/sum(pg0); %Participation factor

%Update AGC parameters
case_sys = case_sys.updateAGCparameters(Kg,Kpart);

%Save the updated object
filename = 'PNobj.mat';
save(filename,'case_sys');

%% Create Jacobian matrices to improve computational efficiency
alwaysCreate = 0;
[isJacobianMade] = computeJacobianMatrix(case_sys,alwaysCreate);

%% Create the Mass matrix
Mass = blkdiag(eye(4*case_sys.nG),zeros(2*case_sys.nG,2*case_sys.nG),zeros(2*case_sys.nN,2*case_sys.nN));

%Save the mass matrix
filename = 'MassMatrix.mat';
save(filename,'Mass');

%% Create the Mass matrix for AGC
MassAGC = blkdiag(eye(4*case_sys.nG + 1),zeros(2*case_sys.nG,2*case_sys.nG),zeros(2*case_sys.nN,2*case_sys.nN));

%Save the mass matrix
filename = 'MassMatrixAGC.mat';
save(filename,'MassAGC');

%% Simulate the steady-state conditions of the power network

%ODE solver settings
options1 = odeset('Mass',Mass,'MassSingular','yes','MStateDependence','none', ...
            'RelTol',1e-7,'AbsTol',1e-6,'Stats','on');

%Include Jacobian        
options1.Jacobian = @JacobianFunc;

%Initial conditions
x_pre0 = [case_sys.x0; case_sys.a0; case_sys.V0];

%Solve the nonlinear DAE without noise
if isSimulSS == 0
    [t1,x_pre] = ode15s(@(t1,x_pre)powerNetwork_NDAE(t1,x_pre,case_sys.u0,pd0,qd0,case_sys),tspan1,x_pre0,options1);
end

%% Create step disturbance
pl0plus = pl0;
ql0plus = ql0;
pr0plus = pr0;
qr0plus = qr0;
pr0plus(ren_set) = pr0(ren_set) + RenDisturbance.*pl0(ren_set);
pl0plus = pl0plus + LoadDisturbance.*pl0;
pd0plus = pl0plus - pr0plus;
qd0plus = ql0plus - qr0plus;

%% Initialization after disturbance
delta0plus = delta0;
E0plus = E0;
TM0plus = TM0;

%% Solve for initial conditions for generator power and bus voltage
[pg0plus,qg0plus,v0plus,theta0plus] = getInitCondAfterDisturbance(case_sys,delta0plus,E0plus,pd0plus,qd0plus);

%% Update again initial conditions and disturbances to object
case_sys = case_sys.updateInitCondAfterDist(delta0plus,E0plus,TM0plus,pd0plus,qd0plus,pl0plus,ql0plus,pr0plus,qr0plus,...
            pg0plus,qg0plus,v0plus,theta0plus);

%Save the updated object
filename = 'PNobj.mat';
save(filename,'case_sys');

%% Simulate the power network after step disturbance due to power mismatch with ode15s

%Setting up initial conditions for simulation
x0plus = [delta0plus; w0*ones(case_sys.nG,1); E0plus; TM0plus];
a0plus = [pg0plus; qg0plus];
V0plus = [v0plus; theta0plus];  

%Initial conditions
x_preD0 = [x0plus; a0plus; V0plus];

%Compute initial conditions for derivative
[initdot] = case_sys.f_gen_dyn_compact(case_sys.delta0plus,case_sys.omega0plus,case_sys.E0plus,case_sys.TM0plus,...
        u0(1:2:2*case_sys.nG),u0(2:2:2*case_sys.nG),case_sys.pg0plus,case_sys.v0plus,case_sys.theta0plus);
xdot0plus0 = zeros(size(x_preD0));
xdot0plus0(1:4*case_sys.nG,1) = initdot;

%ODE solver settings
options2 = odeset('Mass',Mass,'MassSingular','yes','MStateDependence','none', ...
            'RelTol',1e-7,'AbsTol',1e-6,'Stats','on','InitialSlope',xdot0plus0);

%Include Jacobian        
options2.Jacobian = @JacobianFunc;

%Solve the nonlinear DAE without noise
if isSimulDist15s == 1
    [t2,x_pre2] = ode15s(@(t2,x_pre2)powerNetwork_NDAE(t2,x_pre2,u0,pd0plus,qd0plus,case_sys),tspan1,x_preD0,options2);
    plot(t2,x_pre2(:,4))
end

%% Simulate the power network after step disturbance due to power mismatch with ode15i

%ODE solver settings
options3decic = odeset('Jacobian',@JacobianFunc15i,'RelTol',1e-7,'AbsTol',1e-6,'Stats','on');

%Compute initial conditions
[x_pre20,x_pre2p0] = decic(@(t,x_res,x_resp)f_NDAE_implicit(t,x_res,x_resp,u0,pd0plus,qd0plus,Mass,case_sys)...
    ,0,x_preD0,[],xdot0plus0,[],options3decic);

%ODE solver settings
options3 = odeset('Jacobian',@JacobianFunc15i,'RelTol',1e-7,'AbsTol',1e-6,'Stats','on');

%Solve the nonlinear DAE without noise
if isRunDist15i == 1
[t2,x_pre2] = ode15i(@(t4,x_pre2,x_pre2dot)f_NDAE_implicit(t4,x_pre2,x_pre2dot,u0,pd0plus,qd0plus,Mass,case_sys),...
            tspan1,x_pre20,x_pre2p0,options3);       
plot(t2,x_pre2(:,4))
end

%% Get matrices for nonlinear DAE 
sys.E = blkdiag(case_sys.E_tilD(),zeros(2*case_sys.nG+2*case_sys.nN,2*case_sys.nG+2*case_sys.nN));
sys.A = blkdiag(case_sys.A_tilD(),case_sys.A_tilA());
sys.G = blkdiag(case_sys.G_tilD(),case_sys.G_tilA(),eye(2*case_sys.nN));
sys.Bu = [case_sys.B_tilU(); zeros(2*case_sys.nG+2*case_sys.nN,2*case_sys.nG)];
sys.Bq = [zeros(2*case_sys.nG,2*case_sys.nR+2*case_sys.nLN); case_sys.B_tilQ()];
sys.F = [case_sys.F_tilS(); zeros(2*case_sys.nN,1)];

%Smaller matrices
sys.Ed = case_sys.E_tilD();
sys.Ad = case_sys.A_tilD();
sys.Aa = case_sys.A_tilA();
sys.Gd = case_sys.G_tilD();
sys.Ga = blkdiag(case_sys.G_tilA(),eye(2*case_sys.nN));
sys.Bud = case_sys.B_tilU();

%Get matrices for linear DAE
[Ad_d,Ad_a,Aa_d,Aa_a,Bd,Bq] = getLinearDAEmat(case_sys);

%Get matrices for linear ODE
sys.Alin = full(Ad_d-Ad_a*inv(Aa_a)*Aa_d);
sys.Bulin = full(Bd);
sys.Bqlin = sys.Bq;

%Matrices size
sys.nx = size(sys.A,1);
sys.nu = size(sys.Bu,2);
sys.nq = size(sys.Bq,2);
sys.ng = size(sys.G,2);
sys.re = rank(sys.E);
sys.nxd = size(sys.Ad,1);
sys.nxa = size(sys.Aa,1);
sys.ngd = size(sys.Gd,2);
sys.nga = size(sys.Ga,2);
sys.nud = size(sys.Bud,2);
sys.red = rank(sys.Ed);
sys.nxlin = size(sys.Alin,1);
sys.nulin = size(sys.Bulin,2);

%Bounds on the nonlinearity
sys.H = 2.0*eye(sys.nx);
sys.Hd = sqrt(2*eye(sys.nxd));
sys.Ha = sqrt(2*eye(sys.nxa));

%% Compute a stabilizing controller gain matrix
%Nonlinear DAE control
%[K] = NDAE_control(sys); Kd = K(:,1:4*case_sys.nG); % becasue we are designing controller for xd only.
% [Kd,K] = NDAE_reduced_control(sys);
% [Kd,K] = NDAE_reduced_control_2(sys);
[Kd,K,ndae_lmi_time] = NDAE_reduced_control_minimized(sys);
% [Kd,K] = NDAE_reduced_control_input_constr(sys,x0plus);

%LQR design
[Kdlqr,Klqr,lqr_lmi_time] = ODE_LQR_reduced_control(sys,x0plus);

%Construct K for AGC with the aid from NDAE control/LQR for the exciter dynamics
Kdagc = Kdlqr;
Kdagc(case_sys.nG+1:end,1:end) = zeros(case_sys.nG,size(Kdagc,2));
temp = [zeros(case_sys.nG,1); Kpart];
Kdagc = [Kdagc temp];
Kagc = [Kdagc zeros(sys.nud,sys.nxa)];

%Update controller
case_sys = case_sys.updateControlGain(Kd,K,Kdlqr,Klqr,Kdagc,Kagc);

%Save the updated object
filename = 'PNobj.mat';
save(filename,'case_sys');

%% Simulate the power network after step disturbance with controller

%ODE solver settings
options4 = odeset('Mass',sys.E,'MassSingular','yes','MStateDependence','none', ...
            'RelTol',1e-6,'AbsTol',1e-6,'Stats','on');

%Include Jacobian        
options4.Jacobian = @JacobianFuncControl;

%Compute initial inputs after disturbance
u0plus = controlInput(case_sys,sys,u0,x_preD0,Kd);

%Compute initial conditions for derivative
[initdot] = case_sys.f_gen_dyn_compact(case_sys.delta0plus,case_sys.omega0plus,case_sys.E0plus,case_sys.TM0plus,...
        u0plus(1:2:2*case_sys.nG),u0plus(2:2:2*case_sys.nG),case_sys.pg0plus,case_sys.v0plus,case_sys.theta0plus);
xdot0plus = zeros(size(x_preD0));
xdot0plus(1:4*case_sys.nG,1) = initdot;
options4.InitialSlope = xdot0plus;

%Solve the nonlinear DAE with control without noise
if isRunDist15s == 1
    [t3,x_res] = ode15s(@(t3,x_res)powerNetwork_NDAE(t3,x_res,controlInput(case_sys,sys,u0,x_res,Kd),...
                    pd0plus,qd0plus,case_sys),tspan1,x_preD0,options4);
stateidx = 4;
h(1) = figure;
plot(t2,x_pre2(:,stateidx),'r',t3,x_res(:,stateidx),'b--');
h(2) = figure ;
plot(t2,x_pre2(:,stateidx+1),'r',t3,x_res(:,stateidx+1),'b--');
h(3) = figure ;
plot(t2,x_pre2(:,stateidx+2),'r',t3,x_res(:,stateidx+2),'b--');
end

%% Simulate the power network after step disturbance with controller with ode15i

%ODE solver settings
options5decic = odeset('Jacobian',@JacobianFuncControl15i,'RelTol',1e-7,'AbsTol',1e-6,'Stats','on');

%Compute initial conditions
[x_out0,x_outp0] = decic(@(t,x_res,x_resp)f_NDAE_implicit(t,x_res,x_resp,controlInput(case_sys,sys,u0,x_preD0,Kd),pd0plus,qd0plus,Mass,case_sys)...
    ,0,x_preD0,[],xdot0plus,[],options5decic);

%ODE solver settings
options5 = odeset('Jacobian',@JacobianFuncControl15i,'RelTol',1e-7,'AbsTol',1e-6,'Stats','on');

%Solve the ODE
if isRunDist15i == 1
[t4,x_res2] = ode15i(@(t4,x_res2,x_resdot)f_NDAE_implicit(t4,x_res2,x_resdot,controlInput(case_sys,sys,u0,x_res2,Kd),pd0plus,qd0plus,Mass,case_sys),...
            tspan1,x_out0,x_outp0,options5);
end

%% Simulate the LQR-based power network after step disturbance with controller with ode15i

%Compute initial inputs after disturbance
u0plusLQR = controlInput(case_sys,sys,u0,x_preD0,Kdlqr);

%Compute initial conditions for derivative
[initdot] = case_sys.f_gen_dyn_compact(case_sys.delta0plus,case_sys.omega0plus,case_sys.E0plus,case_sys.TM0plus,...
        u0plusLQR(1:2:2*case_sys.nG),u0plusLQR(2:2:2*case_sys.nG),case_sys.pg0plus,case_sys.v0plus,case_sys.theta0plus);
xdot0plusLQR = zeros(size(x_preD0));
xdot0plusLQR(1:4*case_sys.nG,1) = initdot;

%ODE solver settings
options5decic = odeset('Jacobian',@JacobianFuncControl15i,'RelTol',1e-7,'AbsTol',1e-6,'Stats','on');

%Compute initial conditions
[x_out0,x_outp0] = decic(@(t,x_res,x_resp)f_NDAE_implicit(t,x_res,x_resp,controlInput(case_sys,sys,u0,x_preD0,Kdlqr),pd0plus,qd0plus,Mass,case_sys)...
    ,0,x_preD0,[],xdot0plusLQR,[],options5decic);

%ODE solver settings
options5 = odeset('Jacobian',@JacobianFuncControl15i,'RelTol',1e-7,'AbsTol',1e-6,'Stats','on');

%Solve the ODE
if isRunDistLQR15i == 1
[t5,x_res3] = ode15i(@(t5,x_res3,x_resdot)f_NDAE_implicit(t5,x_res3,x_resdot,controlInput(case_sys,sys,u0,x_res3,Kdlqr),pd0plus,qd0plus,Mass,case_sys),...
            tspan1,x_out0,x_outp0,options5);
end

%% Simulate the AGC-LQR-based power network after step disturbance with controller with ode15i

%Compute y0
y0 = 0; %ones(1,case_sys.nG)*(case_sys.Pg0plus-case_sys.Pg0);

%Setting up initial conditions for simulation
x0plusAGC = [delta0plus; w0*ones(case_sys.nG,1); E0plus; TM0plus; y0];

%Initial conditions
x_preD0AGC = [x0plusAGC; a0plus; V0plus];

%Compute initial inputs after disturbance
u0plusAGC = controlInputAGC(case_sys,sys,u0,x_preD0AGC,Kdagc,x0plusAGC);

%Compute initial conditions for derivative
[initdot] = case_sys.f_gen_dyn_AGC(case_sys.delta0plus,case_sys.omega0plus,case_sys.E0plus,case_sys.TM0plus,...
        u0plusAGC(1:2:2*case_sys.nG),u0plusAGC(2:2:2*case_sys.nG),case_sys.pg0plus,case_sys.v0plus,case_sys.theta0plus,y0);
xdot0plusAGC = zeros(size(x_preD0AGC));
xdot0plusAGC(1:4*case_sys.nG + 1,1) = initdot;

%ODE solver settings
% options6decic = odeset('Jacobian',@JacobianFuncAGC15i,'RelTol',1e-7,'AbsTol',1e-6,'Stats','on');
options6decic = odeset('RelTol',1e-7,'AbsTol',1e-6,'Stats','on');

%Compute initial conditions
[x_out0,x_outp0] = decic(@(t,x_res,x_resp)f_NDAE_AGC_implicit(t,x_res,x_resp,controlInputAGC(case_sys,sys,u0,x_preD0AGC,Kdagc,x0plusAGC),...
    pd0plus,qd0plus,MassAGC,case_sys),0,x_preD0AGC,[],xdot0plusAGC,[],options6decic);

%ODE solver settings
% options6 = odeset('Jacobian',@JacobianFuncAGC15i,'RelTol',1e-7,'AbsTol',1e-6,'Stats','on');
options6 = odeset('RelTol',1e-7,'AbsTol',1e-6,'Stats','on');

%Solve the ODE
if isRunDistLQR15i == 1
[t6,x_res4] = ode15i(@(t6,x_res4,x_resdot)f_NDAE_AGC_implicit(t6,x_res4,x_resdot,controlInputAGC(case_sys,sys,u0,x_res4,Kdagc,x0plusAGC),...
    pd0plus,qd0plus,MassAGC,case_sys),tspan1,x_out0,x_outp0,options6);
end

stateidx = 4;
h(7) = figure;
plot(t2,x_pre2(:,stateidx),'r',t4,x_res2(:,stateidx),'b--',t5,x_res3(:,stateidx),'g-.',t6,x_res4(:,stateidx),'k:');
legend('no control','NDAE','LQR','AGC');
set(gcf,'color','w');
h(8) = figure ;
plot(t2,x_pre2(:,stateidx+1),'r',t4,x_res2(:,stateidx+1),'b--',t5,x_res3(:,stateidx+1),'g-.',t6,x_res4(:,stateidx+1),'k:');
legend('no control','NDAE','LQR','AGC');
set(gcf,'color','w');
h(9) = figure ;
plot(t2,x_pre2(:,stateidx+2),'r',t4,x_res2(:,stateidx+2),'b--',t5,x_res3(:,stateidx+2),'g-.',t6,x_res4(:,stateidx+2),'k:');
legend('no control','NDAE','LQR','AGC');
set(gcf,'color','w');

% %% Save Simulation Data
% file_name_data = sprintf('Data_scenario2_%s_dist_%s.mat',casenameraw,distChar);
% save(file_name_data,'casenameraw','case_sys','sys','t2','x_pre2','t4','x_res2','t5','x_res3','t6','x_res4',...
%     'Kd','K','Kdlqr','Klqr','Kdagc','Kagc','ndae_lmi_time','lqr_lmi_time','distChar','RenDisturbance','LoadDisturbance');

%% Collection of Functions

function fx = f_NDAE_AGC_implicit(t,x,x_dot,u,pd,qd,Mass,case_sys)
    %Dynamics
    fx = -Mass*x_dot + powerNetwork_NDAE_AGC(t,x,u,pd,qd,case_sys);
end

function fx = f_NDAE_implicit(t,x,x_dot,u,pd,qd,Mass,case_sys)
    %Dynamics
    fx = -Mass*x_dot + powerNetwork_NDAE(t,x,u,pd,qd,case_sys);
end

function [phaseOut] = degrees2radians(phaseIn)
%DEGREES2RADIANS converts degrees to radians
%   degrees2radians( phaseIn ) converts the angle phaseIn in degrees to its
%   equivalent in radians.
% Created by Hafez Bazrafshan 8/26/2016
phaseOut=-pi+(pi/180)*(phaseIn+180);
end
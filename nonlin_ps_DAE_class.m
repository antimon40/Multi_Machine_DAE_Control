classdef nonlin_ps_DAE_class
% Class for power systems DAE 
% Author: Sebastian A. Nugroho
% Date: 1/28/2021

   properties
       %Network parameters
       casename; %Power system case name
       nN; %Number of buses (N = G + L)
       nG; %Number of generator buses
       nL; %Number of load buses
       nM; %Number of buses with PMU measurements
       Ybus; %Bus admittance matrix
       Gbus;
       Bbus;
       Yf; %From branch admittance matrix
       Yt; %To branch admittance matrix
       Yft; %branch admittance matrix
       node_set; %The set of node labels for the network and has size(nN,1)
       gen_set; %The set of generator labels in the network and has size(nG,1)
       load_set; %The set of load labels with no generators and has size(nL,1) (may include zero loads)
       Sbase; %Base MVA
       w0; %Synchronous frequency
       Cg; %Generators indicator matrix of size(N,G)
       Cr; %Renewables indicator matrix of size(N,R)
       Cl; %Loads indicator matrix of size(N,L)
       ren_set; %The set of renewables
       load_set_nonz; %The set of NONZERO loads
       nR; %Number of renewables
       nLN; %Number of NONZERO loads
       
       Genset; %Set of generators data
       
       %Generator parameters
       M;
       D;
       xq;
       xd;
       xqp;
       xdp;
       Tq0p;
       Td0p;
       Rs;
       Tch;
       Rd;
       
       %Initial conditions
       v0;
       theta0;
       pg0;
       qg0;
       pl0;
       ql0;
       pd0;
       qd0;
       pr0;
       qr0;

       %Set mapping
       GenSetMap;
       LoadSetMap;
       
       %PF solutions 
       isPFok;
       
       %Generator dynamic IC
       delta0;
       Ep0;
       TM0;
       EFd0;
       Tr0;
       
       %DAE initial conditions (at steady-state)
       x0;
       a0;
       V0;
       u0;
       
       %Initial conditions after disturbance
       pd0plus;
       qd0plus;
       delta0plus;
       omega0plus;
       E0plus;
       TM0plus;
       pg0plus;
       qg0plus;
       v0plus;
       theta0plus;
       
       %Renewables and loads after disturbance
       pl0plus;
       ql0plus;
       pr0plus;
       qr0plus;
       
       %Stabilizing controller gain
       K;
       Ksparse;
       Kd;
       Kdfull;
       Klqr;
       Klqrfull;
       Kagc;
       Kagcfull;
       
       %AGC constants
       Kg; %Integrator gain
       Kpart; %Vector of participation factors
       
   end
   methods
       %Object constructor
        function obj = nonlin_ps_DAE_class(casename,nN,nG,nL,nM,Ybus,Gbus,Bbus,Yf,Yt,Cg,node_set,M,D,xq,xd,xqp,xdp,Tq0p,Td0p,Rs,Tch,Rd,...
                gen_set,load_set,Sbase,w0,Genset,v0,theta0,pg0,qg0,pl0,ql0,pd0,qd0,pr0,qr0,GenSetMap,LoadSetMap,isPFok,delta0,E0,TM0,EFd0,Tr0,x0,u0,a0,V0,...
                Cr,Cl,ren_set,load_set_nonz,nR,nLN)
            
            obj.casename = casename;
            obj.nN = nN;
            obj.nG = nG;
            obj.nL = nL;
            obj.nM = nM;
            obj.Ybus = Ybus;
            obj.Gbus = Gbus;
            obj.Bbus = Bbus;
            obj.Yf = Yf;
            obj.Yt = Yt;
            obj.Yft = [Yf; Yt];
            obj.Cg = Cg;
            obj.node_set = node_set;
            obj.gen_set = gen_set;
            obj.load_set = load_set;
            obj.Sbase = Sbase;
            obj.w0 = w0;
            obj.Genset = Genset;
            obj.v0 = v0;
            obj.theta0 = theta0;
            obj.pg0 = pg0;
            obj.qg0 = qg0;
            obj.pl0 = pl0;
            obj.ql0 = ql0;
            obj.pd0 = pd0;
            obj.qd0 = qd0;
            obj.pr0 = pr0;
            obj.qr0 = qr0;
            obj.GenSetMap = GenSetMap;
            obj.LoadSetMap = LoadSetMap;
            obj.isPFok = isPFok;
            obj.delta0 = delta0;
            obj.Ep0 = E0;
            obj.TM0 = TM0;
            obj.EFd0 = EFd0;
            obj.Tr0 = Tr0;
            obj.x0 = x0;
            obj.a0 = a0;
            obj.V0 = V0;
            obj.u0 = u0;
            obj.M  = M;
            obj.D = D;
            obj.xq =xq;
            obj.xd = xd;
            obj.xqp = xqp;
            obj.xdp = xdp;
            obj.Tq0p =Tq0p;
            obj.Td0p = Td0p;
            obj.Rs = Rs;
            obj.Tch = Tch;
            obj.Rd = Rd;
            obj.Cr = Cr; 
            obj.Cl = Cl; 
            obj.ren_set = ren_set; 
            obj.load_set_nonz = load_set_nonz; 
            obj.nR = nR; 
            obj.nLN = nLN; 

        end
        
        %Update initial conditions after disturbance
        function obj = updateInitCondAfterDist(obj,delta0plus,E0plus,TM0plus,pd0plus,qd0plus,pl0plus,ql0plus,pr0plus,qr0plus,...
                pg0plus,qg0plus,v0plus,theta0plus)
            
            obj.pd0plus = pd0plus;
            obj.qd0plus = qd0plus;
            obj.delta0plus = delta0plus;
            obj.E0plus = E0plus;
            obj.TM0plus = TM0plus;
            obj.pl0plus = pl0plus;
            obj.ql0plus = ql0plus;
            obj.pr0plus = pr0plus;
            obj.qr0plus = qr0plus;
            obj.pg0plus = pg0plus;
            obj.qg0plus = qg0plus;
            obj.v0plus = v0plus;
            obj.theta0plus = theta0plus;
            obj.omega0plus = obj.w0*ones(obj.nG,1);
        end
        
        %Update controller gain matrix
        function obj = updateControlGain(obj,Kd,Kdfull,Klqr,Klqrfull,Kagc,Kagcfull)
            obj.Kd = Kd;
            obj.Kdfull = Kdfull;
            obj.Klqr = Klqr;
            obj.Klqrfull = Klqrfull;
            obj.Kagc = Kagc;
            obj.Kagcfull = Kagcfull;
        end
        
        %Update AGC constants
        function obj = updateAGCparameters(obj,Kg,Kpart)
            obj.Kg = Kg; 
            obj.Kpart = Kpart; 
        end
        
        %Generator nonlinear dynamics model (non state-space form)
        function [delta_dot,omega_dot,Ep_dot,Tm_dot] = f_gen_dyn(obj,delta,omega,Ep,Tm,Efd,Tr,Pg,v,theta)
            %Inputs: generator's states delta, omega, Ep, Tm, inputs Efd, Tr, active power Pg, bus voltage v, theta
            %Outputs: time derivative of the states
            
            if isnumeric(delta) && isnumeric(omega) && isnumeric(Ep) && isnumeric(Tm)
                delta_dot = zeros(obj.nG,1); omega_dot = zeros(obj.nG,1); Ep_dot = zeros(obj.nG,1); Tm_dot = zeros(obj.nG,1); 
            else
                delta_dot = sym(zeros(obj.nG,1)); omega_dot = sym(zeros(obj.nG,1)); Ep_dot = sym(zeros(obj.nG,1)); Tm_dot = sym(zeros(obj.nG,1));
            end
            
            %Generator's differential equations
            delta_dot(1:end) = omega - obj.w0*ones(obj.nG,1);
            omega_dot(1:end) = (diag(obj.M))\(Tm - Pg - diag(obj.D)*(omega - obj.w0*ones(obj.nG,1)));
            Ep_dot(1:end) = (diag(obj.Td0p))\(-diag((obj.xd./obj.xdp))*Ep + Efd ... 
                    + diag((obj.xd-obj.xdp)./obj.xdp)*(v(obj.gen_set).*cos(delta-theta(obj.gen_set))));
            Tm_dot(1:end) = (diag(obj.Tch))\(-Tm + Tr - diag(1./obj.Rd)*(omega - obj.w0*ones(obj.nG,1)));
        end
        
        %Generator nonlinear dynamics model (non state-space COMPACT form)
        function [out] = f_gen_dyn_compact(obj,delta,omega,Ep,Tm,Efd,Tr,Pg,v,theta)
            [delta_dot,omega_dot,Ep_dot,Tm_dot] = obj.f_gen_dyn(delta,omega,Ep,Tm,Efd,Tr,Pg,v,theta);
            out = [delta_dot; omega_dot; Ep_dot; Tm_dot];
        end

        %Generator's power equations (non state-space form)
        function [val1,val2] = f_gen_alg(obj,delta,Ep,Pg,Qg,v,theta)
            %Inputs: generator's states delta, Ep, active power Pg, reactive power Qg, bus voltage v, theta
            %Outputs: vectors of 2*nG which should be zero if satisfied
            
            if isnumeric(delta) && isnumeric(Ep) && isnumeric(Pg) && isnumeric(Qg)
                val1 = zeros(obj.nG,1); val2 = zeros(obj.nG,1);
            else
                val1 = sym(zeros(obj.nG,1)); val2 = sym(zeros(obj.nG,1)); 
            end
            
            %Generator's algebraic equations
            val1(1:end) = -Pg + diag(1./obj.xdp)*(Ep.*v(obj.gen_set).*sin(delta-theta(obj.gen_set))) ...
                   - diag((obj.xq-obj.xdp)./(2*(obj.xq.*obj.xdp)))*(v(obj.gen_set).^2.*sin(2*(delta-theta(obj.gen_set))));
            val2(1:end) = -Qg + diag(1./obj.xdp)*(Ep.*v(obj.gen_set).*cos(delta-theta(obj.gen_set))) ...
                   - diag((obj.xq+obj.xdp)./(2*(obj.xq.*obj.xdp)))*(v(obj.gen_set).^2) ...
                   - diag((obj.xq-obj.xdp)./(2*(obj.xq.*obj.xdp)))*(v(obj.gen_set).^2.*cos(2*(delta-theta(obj.gen_set))));
        end
        
        %Generator's power equations (non state-space COMPACT form)
        function [out] = f_gen_alg_compact(obj,delta,Ep,Pg,Qg,v,theta)
            [val1,val2] = obj.f_gen_alg(delta,Ep,Pg,Qg,v,theta);
            out = [val1; val2];
        end
        
        %Network's power balance equations (non state-space form)  without loads nor renewables
        function [valGp,valGq,valLp,valLq] = f_pf(obj,Pg,Qg,v,theta)
            %Inputs: active power Pg, reactive power Qg, bus voltage v, theta
            %Outputs: vectors of 2*nN which should be zero if satisfied
            
            if isnumeric(v) && isnumeric(theta) && isnumeric(Pg) && isnumeric(Qg)
                valGp = zeros(obj.nG,1); valGq = zeros(obj.nG,1); valLp = zeros(obj.nL,1); valLq = zeros(obj.nL,1);
            else
                valGp = sym(zeros(obj.nG,1)); valGq = sym(zeros(obj.nG,1)); valLp = sym(zeros(obj.nL,1)); valLq = sym(zeros(obj.nL,1));
            end
            
            %Pre-compute variables for efficiency
            Dv = diag(v);
            Dctheta = diag(cos(theta));
            Dstheta = diag(sin(theta));
            vctheta = cos(theta);
            vstheta = sin(theta);
            
            %Power balance equations
            vala = -obj.Cg*Pg + Dv*Dctheta*obj.Gbus*Dv*vctheta + Dv*Dstheta*obj.Gbus*Dv*vstheta ...
                   + Dv*Dstheta*obj.Bbus*Dv*vctheta - Dv*Dctheta*obj.Bbus*Dv*vstheta;
            valb = -obj.Cg*Qg + Dv*Dstheta*obj.Gbus*Dv*vctheta - Dv*Dctheta*obj.Gbus*Dv*vstheta ...
                   - Dv*Dstheta*obj.Bbus*Dv*vstheta - Dv*Dctheta*obj.Bbus*Dv*vctheta;
               
            %Output
            valGp(1:end) = vala(obj.gen_set);
            valGq(1:end) = valb(obj.gen_set);
            valLp(1:end) = vala(obj.load_set);
            valLq(1:end) = valb(obj.load_set);
        end
        
        %Network's power balance equations (non state-space form)  WITH loads or renewables
        function [valGp,valGq,valLp,valLq] = f_pf_loads(obj,Pg,Qg,v,theta,Pl,Ql)
            %Inputs: active power Pg, reactive power Qg, bus voltage v,
            %theta, loads power Pl, Ql
            %Outputs: vectors of 2*nN which should be zero if satisfied
            
            if isnumeric(v) && isnumeric(theta) && isnumeric(Pg) && isnumeric(Qg)
                valGp = zeros(obj.nG,1); valGq = zeros(obj.nG,1); valLp = zeros(obj.nL,1); valLq = zeros(obj.nL,1);
            else
                valGp = sym(zeros(obj.nG,1)); valGq = sym(zeros(obj.nG,1)); valLp = sym(zeros(obj.nL,1)); valLq = sym(zeros(obj.nL,1));
            end
            
            %Pre-compute variables for efficiency
            Dv = diag(v);
            Dctheta = diag(cos(theta));
            Dstheta = diag(sin(theta));
            vctheta = cos(theta);
            vstheta = sin(theta);
            
            %Power balance equations
            vala = Pl - obj.Cg*Pg + Dv*Dctheta*obj.Gbus*Dv*vctheta + Dv*Dstheta*obj.Gbus*Dv*vstheta ...
                   + Dv*Dstheta*obj.Bbus*Dv*vctheta - Dv*Dctheta*obj.Bbus*Dv*vstheta;
            valb = Ql - obj.Cg*Qg + Dv*Dstheta*obj.Gbus*Dv*vctheta - Dv*Dctheta*obj.Gbus*Dv*vstheta ...
                   - Dv*Dstheta*obj.Bbus*Dv*vstheta - Dv*Dctheta*obj.Bbus*Dv*vctheta;
               
            %Output
            valGp(1:end) = vala(obj.gen_set);
            valGq(1:end) = valb(obj.gen_set);
            valLp(1:end) = vala(obj.load_set);
            valLq(1:end) = valb(obj.load_set);
        end
        
        %Network's power balance equations (non state-space COMPACT form)  WITH loads or renewables
        function [out] = f_pf_loads_compact(obj,Pg,Qg,v,theta,Pl,Ql)
            [valGp,valGq,valLp,valLq] = obj.f_pf_loads(Pg,Qg,v,theta,Pl,Ql);
            out = [valGp; valGq; valLp; valLq];
        end
        
        %Generator nonlinear dynamics model (non state-space form) with AGC
        %dynamics (after disturbance)
        function [out] = f_gen_dyn_AGC(obj,delta,omega,Ep,Tm,Efd,Tr,Pg,v,theta,y)
            [res] = obj.f_gen_dyn_compact(delta,omega,Ep,Tm,Efd,Tr,Pg,v,theta);
            y_dot = obj.Kg*(-y - (1/obj.nG)*((1./obj.Rd)+obj.D)'*(omega-obj.w0*ones(obj.nG,1)) ...
                    + ones(1,obj.nG)*(Pg-obj.pg0));
            out = [res; y_dot];
        end
        
        %Compute steady state matrices
        %E_tildeD
        function E_tilD_m = E_tilD(obj)
           E_tilD_m = eye(4*obj.nG);
        end
        
        %A_tildeD
        function A_tilD_m = A_tilD(obj)
           A_tilD_m = zeros(4*obj.nG,4*obj.nG);
           A_tilD_m(1:obj.nG,obj.nG+1:2*obj.nG) = eye(obj.nG);
           A_tilD_m(obj.nG+1:2*obj.nG,obj.nG+1:2*obj.nG) = -diag(obj.D./obj.M);
           A_tilD_m(obj.nG+1:2*obj.nG,3*obj.nG+1:end) = diag(1./obj.M);
           A_tilD_m(2*obj.nG+1:3*obj.nG,2*obj.nG+1:3*obj.nG) = (diag(obj.Td0p))\(-diag((obj.xd./obj.xdp)));
           A_tilD_m(3*obj.nG+1:end,obj.nG+1:2*obj.nG) = -(diag(obj.Tch))\(diag(1./obj.Rd));
           A_tilD_m(3*obj.nG+1:end,3*obj.nG+1:end) = -diag(1./obj.Tch);
        end
        
        %G_tildeD
        function G_tilD_m = G_tilD(obj)
           G_tilD_m = zeros(4*obj.nG,2*obj.nG);
           G_tilD_m(obj.nG+1:2*obj.nG,1:obj.nG) = -diag(1./obj.M);
           G_tilD_m(2*obj.nG+1:3*obj.nG,obj.nG+1:end) = (diag(obj.Td0p))\(diag((obj.xd-obj.xdp)./obj.xdp));
        end
        
        %B_tildeU
        function B_tilU_m = B_tilU(obj)
           B_tilU_m = zeros(4*obj.nG,2*obj.nG);
           B_tilU_m(2*obj.nG+1:3*obj.nG,1:obj.nG) = diag(1./obj.Td0p);
           B_tilU_m(3*obj.nG+1:end,obj.nG+1:end) = diag(1./obj.Tch);
        end
        
        %F_tildeS
        function F_tilS_m = F_tilS(obj)
           F_tilS_m = [-ones(obj.nG,1); obj.D./obj.M; zeros(obj.nG,1); 1./(obj.Rd.*obj.Tch)];
        end
        
        %G_tildeA
        function G_tilA_m = G_tilA(obj)
           G_tilA_m = zeros(2*obj.nG,5*obj.nG);
           G_tilA_m(1:obj.nG,1:obj.nG) = diag(1./obj.xdp);
           G_tilA_m(1:obj.nG,obj.nG+1:2*obj.nG) = -diag((obj.xq-obj.xdp)./(2*(obj.xq.*obj.xdp)));
           G_tilA_m(obj.nG+1:end,2*obj.nG+1:3*obj.nG) = diag(1./obj.xdp);
           G_tilA_m(obj.nG+1:end,3*obj.nG+1:4*obj.nG) = -diag((obj.xq+obj.xdp)./(2*(obj.xq.*obj.xdp)));
           G_tilA_m(obj.nG+1:end,4*obj.nG+1:end) = -diag((obj.xq-obj.xdp)./(2*(obj.xq.*obj.xdp)));
        end
        
        %A_tildeA
        function A_tilA_m = A_tilA(obj)
            A_tilA_m = zeros(2*obj.nG+2*obj.nN,2*obj.nG+2*obj.nN);
            A_tilA_m(1:2*obj.nG,1:2*obj.nG) = -eye(2*obj.nG);
            A_tilA_m(2*obj.nG+1:4*obj.nG,1:2*obj.nG) = -eye(2*obj.nG);
        end
        
        %B_tildeQ
        function B_tilQ_m = B_tilQ(obj)
            B_tilQ_m = [];
            
            %Pr and Qr
            B_tilQGP = zeros(obj.nG,length(obj.ren_set));
            B_tilQGQ = zeros(obj.nG,length(obj.ren_set));
            B_tilQLP = zeros(obj.nN,length(obj.ren_set));
            B_tilQLQ = zeros(obj.nN,length(obj.ren_set));
            for i = 1:length(obj.ren_set)
                if ismember(obj.ren_set(i),obj.gen_set) %Generator
                    idxt = find(obj.gen_set == obj.ren_set(i));
                    B_tilQGP(idxt,i) = -1;
                    B_tilQGQ(idxt,i) = -1;
                elseif ismember(obj.ren_set(i),obj.load_set) %Load
                    idxt = find(obj.load_set == obj.ren_set(i));
                    B_tilQLP(idxt,i) = -1;
                    B_tilQLQ(idxt,i) = -1;
                end
            end
            %Cut matrices for loads
            B_tilQLP(obj.gen_set,:) = [];
            B_tilQLQ(obj.gen_set,:) = [];
            
            %Pl and Ql
            B_tilQGP2 = zeros(obj.nG,length(obj.load_set_nonz));
            B_tilQGQ2 = zeros(obj.nG,length(obj.load_set_nonz));
            B_tilQLP2 = zeros(obj.nN,length(obj.load_set_nonz));
            B_tilQLQ2 = zeros(obj.nN,length(obj.load_set_nonz));
            for i = 1:length(obj.load_set_nonz)
                if ismember(obj.load_set_nonz(i),obj.gen_set) %Generator
                    idxt = find(obj.gen_set == obj.load_set_nonz(i));
                    B_tilQGP2(idxt,i) = 1;
                    B_tilQGQ2(idxt,i) = 1;
                elseif ismember(obj.load_set_nonz(i),obj.load_set) %Load
                    idxt = find(obj.load_set == obj.load_set_nonz(i));
                    B_tilQLP2(idxt,i) = 1;
                    B_tilQLQ2(idxt,i) = 1;
                end
            end
            %Cut matrices for loads
            B_tilQLP2(obj.gen_set,:) = [];
            B_tilQLQ2(obj.gen_set,:) = [];
            
            B_tilQ_m = [blkdiag(B_tilQGP,B_tilQGQ) blkdiag(B_tilQGP2,B_tilQGQ2); ...
                        blkdiag(B_tilQLP,B_tilQLQ) blkdiag(B_tilQLP2,B_tilQLQ2)];
        end

    end
end
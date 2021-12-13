classdef psgen_class
% Class for generator, 4-th order nonlinear DAE model
% By: Sebastian A. Nugroho
% Date: 7/22/2020
   properties
       %Generator's dynamic parameters
       w0; %Synchronous frequency
       M; %Rotor's inertia constant
       D; %Damping coefficient
       xq; %Direct axis synchronous reactance in q-axis 
       xd; %Direct axis synchronous reactance in d-axis 
       xqp; %Direct axis transient reactance in q-axis 
       xdp; %Direct axis transient reactance in d-axis 
       Tq0p; %Open circuit time constant in q-axis
       Td0p; %Open circuit time constant in d-axis
       Tch; %Chest valve time constants
       Rd; %Regulation constants
       
       %Generator's algebraic parameters
       Rs; %Stator resistance
       
       %Generator bus index
       bus_idx;
       
       %Generator index
       gen_idx;
       
       %Generator base MVA
       genSbase;
   end
   methods
        function obj = psgen_class(w0,M,D,xq,xd,xqp,xdp,Tq0p,Td0p,Tch,Rd,Rs,bus_idx,gen_idx,genSbase)
%             if nargin == 13
                obj.w0 = w0;
                obj.M = M;
                obj.D = D;
                obj.xq = xq;
                obj.xd = xd;
                obj.xqp = xqp;
                obj.xdp = xdp;
                obj.Tq0p = Tq0p;
                obj.Td0p = Td0p;
                obj.Tch = Tch;
                obj.Rd = Rd;
                obj.Rs = Rs;
                obj.bus_idx = bus_idx;
                obj.gen_idx = gen_idx;
                obj.genSbase = genSbase;
%             else
%                 disp('Number of arguments does not match. Please check.')
%             end
        end
    end
end
function [w0,M,D,xq,xd,xqp,xdp,Tq0p,Td0p,Rs,Sbase2] = getNetworkDescriptionPST(net_desc)
% Get power system's network description from PST
% input: net_desc is the case file name for PST
% outputs:
% 1. w0 is the synchronous frequency
% 2. M; %Rotor's inertia constant
% 3. D; %Damping coefficient
% 4. xq; %Direct axis synchronous reactance in q-axis 
% 5. xd; %Direct axis synchronous reactance in d-axis 
% 6. xqp; %Direct axis transient reactance in q-axis 
% 7. xdp; %Direct axis transient reactance in d-axis 
% 8. Tq0p; %Open circuit time constant in q-axis
% 9. Td0p; %Open circuit time constant in d-axis
% 10. Rs; %Stator resistance
% 11. Sbase2: base MVA from PST
% Author: Sebastian A. Nugroho
% Date: 7/22/2020

%Load data
lfile =length(net_desc);
% strip off .m and convert to lower case
net_desc = lower(net_desc(1:lfile-2));
eval(net_desc);

%Assign generator data
sys_freq = 60; %60 Hz
w0 = 2*pi*sys_freq; %rad/sec
% w0 = 1; %pu
% mac_pot = 100/mac_con(mac_num,3);
M = 0.01*mac_con(:,16);
%M = mac_con(:,16)/(pi*60);
xd = 1.2*mac_con(:,6);
xdp = mac_con(:,7);
xq = mac_con(:,11);
xqp = mac_con(:,12);
Td0p = mac_con(:,9);
Tq0p = 2*mac_con(:,14);
D = 0.2*mac_con(:,17);
Rs = mac_con(:,5);
Sbase2 = mac_con(:,3);

end
function [nN,nG,nL,Ybus,Gbus,Bbus,node_set,gen_set,load_set,...
    Ysh,Cf,Ct,Cg,Yff,Yft,Ytf,Ytt,Sbase,Yf,Yt] = getNetworkDescriptionMATPOWER(net_desc)
% Get power system's network description from MATPOWER
% Code is developed from 'networkParams.m' by Hafez Bazrafshan
% Input: net_desc is obtained from running loadcase function from MATPOWER
% Outputs:
% 1. nN is the number of nodes (or buses)
% 2. nG is the number of generator buses
% 3. nL is the number of load buses
% 4. Ybus: is the complex bus admittance matrix of size(nN,nN)
% 5. Gbus: is the real part of Ybus. 
% 6. Bbus: is the imaginary part of Ybus. 
% 7. node_set: is the set of node labels for the network and has size(nN,1)
% 8. gen_set: is the set of generator labels in the network and has size(nG,1)
% 9. load_set: is the set of load labels with no generators and has size(nL,1)
% 10. Ysh: is the matrix of node shunt admittances (diagonal) of size(nN,nN)
% 11. Cf: is the from matrix of size(B,N), where B is the number of branches
% 12. Ct: is the to matrix of size(B,N), where B is the number of branches 
% 13. Cg: is the generator indicator matrix of size(nN,nG) 
% 14. Sbase: base MVA 
% 15. Yf: from branch admittance matrix 
% 16. Yt: to branch admittance matrix 
% Author: Sebastian A. Nugroho
% Date: 7/22/2020

%Labels
NodeLabels=union(net_desc.branch(:,1),net_desc.branch(:,2)); 
GenLabels=net_desc.gen(:,1); 
LoadLabels=setdiff(NodeLabels,GenLabels); 

%Number of nodes
nN=length(NodeLabels); 
nG=length(GenLabels); 
nL=length(LoadLabels);

%Set of nodes, generators, and loads
node_set=(1:nN).';
gen_set=zeros(nG,1); 
for ii=1:nG
    gen_set(ii)=find(NodeLabels==GenLabels(ii));
end
load_set=zeros(nL,1);
for ii=1:nL
    load_set(ii)=find(NodeLabels==LoadLabels(ii));
end

%Construct bus and branch admittance matrices
FromNodes=net_desc.branch(:,1); 
ToNodes=net_desc.branch(:,2); 
B=length(FromNodes);  % number of branches

FromNodesIndices=zeros(B,1); 
ToNodesIndices=zeros(B,1);

for ii=1:B
    FromNodesIndices(ii)=find(NodeLabels==FromNodes(ii));
    ToNodesIndices(ii)=find(NodeLabels==ToNodes(ii));
end

Cf=sparse(B,nN); 
Ct=sparse(B,nN); 

Cf(sub2ind([B nN], (1:B).', FromNodesIndices))=1;
Ct(sub2ind([B nN], (1:B).', ToNodesIndices))=1; 

rmn_vec=net_desc.branch(:,3); 
xmn_vec=net_desc.branch(:,4);
zmn_vec=rmn_vec+i*xmn_vec;
ymn_vec=1./zmn_vec;
gmn_vec=real( ymn_vec); 
bmn_vec=imag(ymn_vec);
bcmn_vec=net_desc.branch(:,5); 
taumn_vec=net_desc.branch(:,9); 
taumn_vec(taumn_vec==0)=1;
alphamn_vec=net_desc.branch(:,10);
gs_vec=net_desc.bus(:,5);
bs_vec=net_desc.bus(:,6);
ys_vec=(gs_vec+i*bs_vec)/net_desc.baseMVA; % vector of nodal shunt admittances

YffVec=(gmn_vec+i*(bmn_vec+bcmn_vec./2))./(taumn_vec.^2); 
YftVec=-(gmn_vec+i*bmn_vec)./(taumn_vec.*exp(-i*degrees2radians(alphamn_vec))); 
YtfVec=-(gmn_vec+i*bmn_vec)./(taumn_vec.*exp(i*degrees2radians(alphamn_vec))); 
YttVec=gmn_vec+i*(bmn_vec+bcmn_vec/2); 

Yff=diag(sparse(YffVec)); 
Yft=diag(sparse(YftVec)); 
Ytf=diag(sparse(YtfVec)); 
Ytt=diag(sparse(YttVec)); 
Ysh=diag(sparse(ys_vec)); 

Yf = Yff*Cf + Yft*Ct;
Yt = Ytf*Cf+Ytt*Ct;
Ybus=Cf.'*Yff*Cf+Cf.'*Yft*Ct+Ct.'*Ytf*Cf+Ct.'*Ytt*Ct+Ysh;
Gbus=real(Ybus); 
Bbus=imag(Ybus); 

Cg=sparse(nN,nG); 

Cg(sub2ind([nN nG], gen_set, (1:nG).'))=1;

Sbase = net_desc.baseMVA; 
% omega0 = 2*pi*60; %60 Hz

%Check bus and branch admittance matrices from MATPOWER function
disp('Check MATPOWER data.')
[Ybusm, Yfm, Ytm] = makeYbus(net_desc);

if ~isequal(Ybus,Ybusm)
    disp('Ybus not equal. Use from MATPOWER.')
else
    disp('Ybus OK.')
end

if ~isequal(Yf,Yfm)
    disp('Yf not equal. Use from MATPOWER.')
else
    disp('Yf OK.')
end

if ~isequal(Yt,Ytm)
    disp('Yt not equal. Use from MATPOWER.')
else
    disp('Yt OK.')
end

function [phaseOut] = degrees2radians(phaseIn)
%DEGREES2RADIANS converts degrees to radians
%   degrees2radians( phaseIn ) converts the angle phaseIn in degrees to its
%   equivalent in radians.
% Created by Hafez Bazrafshan 8/26/2016
phaseOut=-pi+(pi/180)*(phaseIn+180);
end

end
function out = SNAIL
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
%------------------------------------------------------------------------------ % t,kmrgd,lam3,S,g1,g2,g3,z1,z2,s1,s2,lam1,lam2,lam4,k1,k2,k3,n1,n2,n3,n4,L0,L1,L2,L3,L4,L5,L6,xnot,nmu,Gam1_0,Gam1_1,Gam1_2,Gam1_3,Gam1_4,Gam1_5,Gam1_6,Gam2_0,Gam2_1,Gam2_2,Gam2_3,Gam2_4,Gam2_5,Gam2_6
function dydt = fun_eval(t,kmrgd,lam3,g1,g2,g3,z1,z2,lam1,k1,k2,k3,n1,n3,L0,L1,L2,L3,L4,L5,L6,xnot,nmu)

%Lm function components 
Mu0=1/(1+kmrgd(1)/xnot)^nmu;
Mu1=(kmrgd(1)/xnot)/(1+kmrgd(1)/xnot)^nmu;
Mu2=(kmrgd(1)/xnot)^2/(1+kmrgd(1)/xnot)^nmu;
Mu3=(kmrgd(1)/xnot)^3/(1+kmrgd(1)/xnot)^nmu;
Mu4=(kmrgd(1)/xnot)^4/(1+kmrgd(1)/xnot)^nmu;
Mu5=(kmrgd(1)/xnot)^5/(1+kmrgd(1)/xnot)^nmu;
Mu6=(kmrgd(1)/xnot)^6/(1+kmrgd(1)/xnot)^nmu;

%Hills functions
H1=(z1^n1+lam1*kmrgd(3)^n1)/(z1^n1+kmrgd(3)^n1);
%H2=(s1^n2+lam2*S^n2)/(s1^n2+S^n2);
H3=(z2^n3+lam3*kmrgd(3)^n3)/(z2^n3+kmrgd(3)^n3);
%H4=(s2^n4+lam4*S^n4)/(s2^n4+S^n4);

%Equations
dydt=[g1*H1 - k1*kmrgd(1);   
    g2*H3 - k2*kmrgd(2);     
    g3*kmrgd(2)*(Mu0*L0 + 6*Mu1*L1 + 15*Mu2*L2 + 20*Mu3*L3 + 15*Mu4*L4 + 6*Mu5*L5 + Mu6*L6) - k3*kmrgd(3)];

%(kmrgd(2)*(Mu0*Gam1_0 + 6*Mu1*Gam1_1 + 15*Mu2*Gam1_2 + 20*Mu3*Gam1_3 + 15*Mu4*Gam1_4 + 6*Mu5*Gam1_5 + Mu6*Gam1_6)) - 
%(kmrgd(2)*(Mu0*Gam2_0 + 6*Mu1*Gam2_1 + 15*Mu2*Gam2_2 + 20*Mu3*Gam2_3 + 15*Mu4*Gam2_4 + 6*Mu5*Gam2_5 + Mu6*Gam2_6)) -
% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(SNAIL);
y0=[0,0,0,0]; 
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% ---------------------------- t,kmrgd,g1,g2,g3,z1,z2,s1,s2,lam1,lam2,lam3,lam4,k1,k2,k3,n1,n2,n3,n4,L0,L1,L2,L3,L4,L5,L6,xnot,nmu,Gam1_0,Gam1_1,Gam1_2,Gam1_3,Gam1_4,Gam1_5,Gam1_6,Gam2_0,Gam2_1,Gam2_2,Gam2_3,Gam2_4,Gam2_5,Gam2_6
function jac = jacobian(t,kmrgd,lam3,g1,g2,g3,z1,z2,lam1,k1,k2,k3,n1,n3,L0,L1,L2,L3,L4,L5,L6,xnot,nmu) %
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,lam3,g1,g2,g3,z1,z2,lam1,k1,k2,k3,n1,n3,L0,L1,L2,L3,L4,L5,L6,xnot,nmu)
%--------------------------------------------------------------------------
function hess = hessians(t,kmrgd,lam3,g1,g2,g3,z1,z2,lam1,k1,k2,k3,n1,n3,L0,L1,L2,L3,L4,L5,L6,xnot,nmu)
%--------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,lam3,g1,g2,g3,z1,z2,lam1,k1,k2,k3,n1,n3,L0,L1,L2,L3,L4,L5,L6,xnot,nmu)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,lam3,g1,g2,g3,z1,z2,lam1,k1,k2,k3,n1,n3,L0,L1,L2,L3,L4,L5,L6,xnot,nmu)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,lam3,g1,g2,g3,z1,z2,lam1,k1,k2,k3,n1,n3,L0,L1,L2,L3,L4,L5,L6,xnot,nmu)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,lam3,g1,g2,g3,z1,z2,lam1,k1,k2,k3,n1,n3,L0,L1,L2,L3,L4,L5,L6,xnot,nmu)
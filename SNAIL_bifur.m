function [x,v,s,h,f] = SNAIL_bifur

curdir = pwd;
init;
cd(curdir)
opt = contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',50000); %50000 %500000
opt=contset(opt,'MinStepsize',1);  %1
opt=contset(opt,'MaxStepsize',100);  %100
opt=contset(opt,'Eigenvalues',1);

% Degradation rate:
k1=0.05; k2=0.5; k3=0.1;       %TI, Both TI&AD

% Transcription rate:
%g1=1200; g2=10; g3=100;      %Both TI&AD
%g1=1500; g2=12.5; g3=100;    %TI only
%g1=1500; g2=15; g3=100;      %Both TI&AD but without SNAIL
g1=1000; g2=65; g3=100;      % Just Toggle Switch


% Hills function threshold :
%z1=120000; z2=24000; s1=175000; s2=180000; xnot=10000;  %Both TI&AD
%z1=200000; z2=50000; s1=180000; s2=180000; xnot=10000;  %TI only
%z1=200000; z2=25000; s1=0; s2=0; xnot=10000;  %Both TI&AD but without SNAIL
z1=200000; z2=50000; %s1=0; s2=0; 
xnot=10000;  %Just Toggle Switch

 % Cooperativity:
n1=2; %n2=0; 
n3=2; %n4=0; 
nmu=6;

% fold change
%lam1=0.1; lam2=0.75; lam3=9.8; lam4=9; %Both TI&AD  %i=3, lam3=9.8 and for i=6, lam3=11.3
%lam1=0.1; lam2=0.4; lam3=15; lam4=10; %TI only 
lam1=0.1; %lam2=0; lam4=0; %Both TI&AD but without SNAIL  
%lam3=70; 
%S=0;

%Translational Inhibition rates:
%L0=1.0; L1=0.6; L2=0.3; L3=0.1; L4=0.1; L5=0.1; L6=0.1;  %TI&AD 
%L0=1.0; L1=0.5; L2=0.2; L3=0.02; L4=0.02; L5=0.02; L6=0;  %L3=0.02; L4=0.02; L5=0.02; L6=0.02 %TI only   
%L0 = 1; L1 = 0.6; L2 = 0.3; L3 = 0.1; L4 = 0.05; L5 = 0.05; L6 = 0.05;    %Both TI&AD but without SNAIL 
L0 = 1; L1 = 0.5; L2 = 0.2; L3 = 0.02; L4 = 0; L5 = 0; L6 = 0; % Just toggle switch

%Active Degradation of miR200 (Gam_1's) and zebMRNA (Gam_2's)
%Gam1_0 = 0; Gam1_1 = 0.007; Gam1_2 = 0.07; Gam1_3 = 0.7; Gam1_4 = 0.7; Gam1_5 = 0.7; Gam1_6 = 0.7;    % Both TI & AD
%Gam2_0 = 0; Gam2_1 = 0.04;  Gam2_2 = 0.2;  Gam2_3 = 1;   Gam2_4 = 1;   Gam2_5 = 1; Gam2_6 = 1;    % Both TI & AD

%%%%%%%%%%% TI only
%Gam1_0 = 0; Gam1_1 = 0;  Gam1_2 = 0;  Gam1_3 = 0;   Gam1_4 = 0;   Gam1_5 = 0; Gam1_6 = 0;    % Both TI & AD
%Gam2_0 = 0; Gam2_1 = 0;  Gam2_2 = 0;  Gam2_3 = 0;   Gam2_4 = 0;   Gam2_5 = 0; Gam2_6 = 0;    % Both TI & AD


%%%%% Both TI & AD WS
%Gam1_0 = 0; Gam1_1 = 0.005; Gam1_2 = 0.05; Gam1_3 = 0.5; Gam1_4 = 0.5; Gam1_5 = 0.5; Gam1_6 = 0.5;     % Both TI & AD WS
%Gam2_0 = 0; Gam2_1 = 0.04;  Gam2_2 = 0.2;  Gam2_3 = 1;   Gam2_4 = 1;   Gam2_5 = 1;   Gam2_6 = 1;       % Both TI & AD WS  

%%%%% Just Toggle Switch
%Gam1_0 = 0; Gam1_1 = 0;  Gam1_2 = 0;  Gam1_3 = 0;   Gam1_4 = 0;   Gam1_5 = 0;   Gam1_6 = 0;     
%Gam2_0 = 0; Gam2_1 = 0;  Gam2_2 = 0;  Gam2_3 = 0;   Gam2_4 = 0;   Gam2_5 = 0;   Gam2_6 = 0;      
                 
% SNAIL   %put S=0 when it is a bifurcation parameter
%S=0;  % S=0

% lam3   %put lam3=0 when it is a bifurcation parameter
lam3=0;

ap = 1; %describes the index of parameter for which the bifurcation is drawn using the init_EP_EP function. Currently, ap=1, thus bifurcation parameter is s (SNAIL levels)
handles = feval(@SNAIL);
tspan = 0:100:5000;   

% initial condition    EP(18721.89 101.15 9892.49) works for all i in both TI&AD 3222 case
            % EP(17837.73925 197.6780136 20267.03457) works for all i in both TI&AD 2121 case, but with lamda 3 changes
x_start = [19952 336 10309];  %[19361 338 150261]; % [2326 2405 1279997]; %[33554.833280 56.5 0];   

%calculating steady state for given initial condition 
[t,x_time] = ode15s(@(t,kmrgd)handles{2}(t,kmrgd,lam3,g1,g2,g3,z1,z2,lam1,k1,k2,k3,n1,n3,L0,L1,L2,L3,L4,L5,L6,xnot,nmu),tspan,x_start); % s1,s2,lam2,lam4,n2,n4,
x_init = [19952; 336; 10309];  %x_time(end,:)'; 

%drawing bifurcation using a continuation method  19743.39 975.8 42074.28
[x0,v0] = init_EP_EP(@SNAIL,x_init,[lam3,g1,g2,g3,z1,z2,lam1,k1,k2,k3,n1,n3,L0,L1,L2,L3,L4,L5,L6,xnot,nmu],ap); % s1,s2,lam2,lam4,,n2,n4,
[x,v,s,h,f] = cont(@equilibrium, x0, v0,opt);


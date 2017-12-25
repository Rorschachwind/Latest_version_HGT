 %{
----------------12/06/2017----------------
This is a record for a former version with 
update part. Not used now.
---------------main function--------------
This method is not finally used, but the update method used in the main
function is useful.
Anyway, this is our improvement and simplification.
This version contain a more complex update part.
%}

clear; clc;
%% global term
global K1 K2 alpha1 alpha2 beta1 beta2 n m Nm mu_G_max mu_R_max mu_Y_max  initConc tspan eta0 A Atype

K1 = 0.1;
K2 = 0.1;
alpha1 = 3.4409e-07;
alpha2 = 0.0017;
beta1 = 25;
beta2 = 25;
n = 2;
m = 2;
Nm = 1E9; % carrying capacity
D=0.15;
xLimit = [0 60];
% antibiotics type 
Atype = 'none';
A = 0; %Antibiotics

if isequal(lower(Atype), 'cm')
        mu_G_max = 0.33;
        mu_R_max = 0.32;
        mu_Y_max = 0.31;
end

if isequal(lower(Atype), 'kan')
        mu_G_max = 0.33;
        mu_R_max = 0.28;
        mu_Y_max = 0.31;
end

tspan= 0.1:0.1:60;
initConc = [1E-3*Nm; 1E-3*Nm; 0];
etaC_prime = 0.0375;
eta0 = (etaC_prime * mu_Y_max / Nm);
%% set up function
% 1:G a R r
% 2: G r R a
% 3: G a R a
% 4: G r R r

%constant etaC without dilution D (the original one)
Version = 0;
[tv,FvS]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,0,Version),tspan,initConc);

%constant etaC with dilution D
Version = 0;
[tv,Fv0]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,Version),tspan,initConc);

% G a R r
Version = 1;
[tv,Fv]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,Version),tspan,initConc);
X =calcX(Fv,Version);
[tv,Fv1]=ode45(@(t,Y) fun_GRY_Hill_D_update(t,Y,D,Version,Fv),tspan, initConc);
[mu_eff_update1,etaGR1,etaYR1]= calcE (Fv1,Version,X);

%G r R a
Version = 2;
[tv,Fv]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,Version),tspan,initConc);
X =calcX(Fv,Version);
[tv,Fv2]=ode45(@(t,Y) fun_GRY_Hill_D_update(t,Y,D,Version,Fv),tspan, initConc);
[mu_eff_update2,etaGR2,etaYR2]= calcE (Fv2,Version,X);

%G a R a
Version = 3;
[tv,Fv]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,Version),tspan,initConc);
X =calcX(Fv,Version);
[tv,Fv3]=ode45(@(t,Y) fun_GRY_Hill_D_update(t,Y,D,Version,Fv),tspan, initConc);
[mu_eff_update3,etaGR3,etaYR3]= calcE (Fv3,Version,X);

%G r R r
Version = 4;
[tv,Fv]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,Version),tspan,initConc);
X =calcX(Fv,Version);
[tv,Fv4]=ode45(@(t,Y) fun_GRY_Hill_D_update(t,Y,D,Version,Fv),tspan, initConc);
[mu_eff_update4,etaGR4,etaYR4]= calcE (Fv4,Version,X);

% eta_total
xx = etaYR1+etaGR1;
yy =  etaYR2+etaGR2;
zz = etaYR3+etaGR3;
ww =  etaYR4+etaGR4;
%% plot All figures
Popfig(tv,FvS,Fv0,Fv1,Fv2,Fv3,Fv4,xx,yy,zz,ww)

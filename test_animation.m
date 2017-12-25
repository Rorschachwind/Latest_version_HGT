%{ 
-------------12/19/2017----------
This is used in this set of function.
This is used to record animation of 
three population relationships.
------------Resistance--------------
         G is sensitive to Cm
         R is resistant to Cm
         Y is resistant to Cm

         G is resistant to Kan
         R is sensitive to Kan
         Y is resistant to Kan
-----Combination of eta_Hill--------
             1: G a R r
             2: G r R a
             3: G a R a
             4: G r R r
%}

%% global term
clear; clc;
global K1 K2 alpha1 alpha2 beta1 beta2 n m Nm mu_G_max mu_R_max mu_Y_max  initConc tspan eta0 A Atype

K1 = 0.01;
K2 = 0.01;
alpha1 = 3.4409e-07;
alpha2 = 0.0017;
beta1 = 25;
beta2 = 25;
n = 2;
m = 2;
Nm = 1E9; % carrying capacity
xLimit = [0 60];
 
%Dilution rate
% D=0.0;

%segregation error (don't change)
kappa=0.0;

% antibiotics type 
Atype = 'none';
A =0.0; 
if isequal(lower(Atype), 'cm')
        mu_G_max = 0.33;
        mu_R_max = 0.32;% cm actually 0.32
        mu_Y_max = 0.31;
elseif isequal(lower(Atype), 'kan')
        mu_G_max = 0.33;
        mu_R_max = 0.28;
        mu_Y_max = 0.31;
elseif isequal(lower(Atype), 'none')
        mu_G_max = 0.33;
        mu_R_max = 0.32;
        mu_Y_max = 0.31;
elseif isequal(lower(Atype), 'both')
        mu_G_max = 0.33;
        mu_R_max = 0.28;
        mu_Y_max = 0.31;    
end

tstep=0.5;
tspan= 0.5:tstep:60;
initConc = [0.3*Nm; 0.3*Nm; 0];
etaC_prime = 0.0375;
eta0 = (etaC_prime * mu_Y_max / Nm);
%% do integral
%{ 
Idea to generate parameter: to keep the eta_C integral as the same as
eta_0, basically can do double integral to the term, or seperately
calculate as followed\
%}

mu_G_eff=0:0.01:mu_G_max;
mu_R_eff=0:0.01:mu_R_max;
mu_Y_eff=0:0.01:mu_Y_max;
Hill_GA = @(mu_G_eff) alpha1 + alpha2 * mu_G_eff.^n./(K1^n+mu_G_eff.^n);
Hill_GR = @(mu_G_eff) alpha1 + alpha2 * K1^n./(K1^n+mu_G_eff.^n);
Hill_RA = @(mu_R_eff) beta1 + beta2 * mu_R_eff.^m./(K2^m+mu_R_eff.^m);
Hill_RR = @(mu_R_eff) beta1 + beta2 * K2^m./(K2^m+mu_R_eff.^m);
Hill_YA=@(mu_Y_eff) alpha1 + alpha2 * mu_Y_eff.^n./(K1^n+mu_Y_eff.^n);
Hill_YR=@(mu_Y_eff) alpha1 + alpha2 * K1.^n./(K1^n+mu_Y_eff.^n);

q_GA=integral(Hill_GA,0,mu_G_max)/mu_G_max;
q_GR=integral(Hill_GR,0,mu_G_max)/mu_G_max;
q_RA=integral(Hill_RA,0,mu_R_max)/mu_R_max;
q_RR=integral(Hill_RR,0,mu_R_max)/mu_R_max;
q_YA=integral(Hill_YA,0,mu_Y_max)/mu_Y_max;
q_YR=integral(Hill_YR,0,mu_Y_max)/mu_Y_max;

%% set up function
%constant etaC without dilution D (the original one)
% Version = 0;
% [tv,FvS]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,0,0,Version,[],[],[]),tspan,initConc);
% DD=[0.01 0.02];
DD=[0.01 0.03 0.05 0.07 0.1 0.12 0.15];
dd=length(DD);
%  AA=[0.1 0.5 1.0 1.5];
% AA=[0.1 0.5 1.0 1.5 2.5 3.0 3.5 4.0 4.5 5.0];
% 2.5 3.0 3.5 4.0 4.5 5.0
% aa=length(AA);
for j =1:dd
    D=DD(j);
    %constant etaC with dilution D
    Version = 0;
    [tv,Fv0]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,[],[],[]),tspan,initConc);
    jj=length(tv);
    Y0_1=Fv0(:,1)/Nm;
    Y0_2=Fv0(:,2)/Nm;
    Y0_3=Fv0(:,3)/Nm;
    % G a R r
    Version = 1;
    [tv,Fv1]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GA,q_RR,q_YA),tspan, initConc);
    [mu_eff_update1,etaGR1,etaYR1]= calcE (Fv1,Version,q_GA,q_RR,q_YA);
    Y1_1=Fv1(:,1)/Nm;
    Y1_2=Fv1(:,2)/Nm;
    Y1_3=Fv1(:,3)/Nm;
    %G r R a
    Version = 2;
    [tv,Fv2]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GR,q_RA,q_YR),tspan, initConc);
    [mu_eff_update2,etaGR2,etaYR2]= calcE (Fv2,Version,q_GR,q_RA,q_YR);
    Y2_1=Fv2(:,1)/Nm;
    Y2_2=Fv2(:,2)/Nm;
    Y2_3=Fv2(:,3)/Nm;
    %G a R a
    Version = 3;
    [tv,Fv3]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GA,q_RA,q_YA),tspan, initConc);
    [mu_eff_update3,etaGR3,etaYR3]= calcE (Fv3,Version,q_GA,q_RA,q_YA);
    Y3_1=Fv3(:,1)/Nm;
    Y3_2=Fv3(:,2)/Nm;
    Y3_3=Fv3(:,3)/Nm;
    %G r R r
    Version = 4;
    [tv,Fv4]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_GR,q_RR,q_YR),tspan, initConc);
    [mu_eff_update4,etaGR4,etaYR4]= calcE (Fv4,Version,q_GR,q_RR,q_YR);
    Y4_1=Fv4(:,1)/Nm;
    Y4_2=Fv4(:,2)/Nm;
    Y4_3=Fv4(:,3)/Nm;

    %% calculate frac and integral
    % constant etaC
    Frac_YR_0=log(1+Y0_3./Y0_2);
    Frac_YR_0=Frac_YR_0(2:jj);
    Intrg_G_0=[];
    Intrg_G_0(1)=(Y0_1(1)+Y0_1(2));
    Intrg_Y_0(1)=(Y0_3(1)+Y0_3(2));
    for i = 2:(jj-1)
        Intrg_G_0(i)=Intrg_G_0(i-1)+Y0_1(i)+Y0_1(i+1);
        Intrg_Y_0(i)=Intrg_Y_0(i-1)+Y0_3(i)+Y0_3(i+1);
    end
    Intrg_G_0=Intrg_G_0.*tstep./2;
    Intrg_Y_0=Intrg_Y_0.*tstep./2;


    % Ga Rr
    Frac_YR_1=log(1+Y1_3./Y1_2);
    Frac_YR_1=Frac_YR_1(2:jj);
    Intrg_G_1=[];
    Intrg_G_1(1)=(Y1_1(1)+Y1_1(2));
    Intrg_Y_1(1)=(Y1_3(1)+Y1_3(2));
    for i = 2:(jj-1)
        Intrg_G_1(i)=Intrg_G_1(i-1)+Y1_1(i)+Y1_1(i+1);
        Intrg_Y_1(i)=Intrg_Y_1(i-1)+Y1_3(i)+Y1_3(i+1);
    end
    Intrg_G_1=Intrg_G_1.*tstep./2;
    Intrg_Y_1=Intrg_Y_1.*tstep./2;

    %Gr Ra
    Frac_YR_2=log(1+Y2_3./Y2_2);
    Frac_YR_2=Frac_YR_2(2:jj);
    Intrg_G_2=[];
    Intrg_G_2(1)=(Y2_1(1)+Y2_1(2));
    Intrg_Y_2(1)=(Y2_3(1)+Y2_3(2));
    for i = 2:(jj-1)
        Intrg_G_2(i)=Intrg_G_2(i-1)+Y2_1(i)+Y2_1(i+1);
        Intrg_Y_2(i)=Intrg_Y_2(i-1)+Y2_3(i)+Y2_3(i+1);
    end
    Intrg_G_2=Intrg_G_2.*tstep./2;
    Intrg_Y_2=Intrg_Y_2.*tstep./2;

    %Ga Ra
    Frac_YR_3=log(1+Y3_3./Y3_2);
    Frac_YR_3=Frac_YR_3(2:jj);
    Intrg_G_3=[];
    Intrg_G_3(1)=(Y3_1(1)+Y3_1(2));
    Intrg_Y_3(1)=(Y3_3(1)+Y3_3(2));
    for i = 2:(jj-1)
        Intrg_G_3(i)=Intrg_G_3(i-1)+Y3_1(i)+Y3_1(i+1);
        Intrg_Y_3(i)=Intrg_Y_3(i-1)+Y3_3(i)+Y3_3(i+1);
    end
    Intrg_G_3=Intrg_G_3.*tstep./2;
    Intrg_Y_3=Intrg_Y_3.*tstep./2;

    %Gr Rr
    Frac_YR_4=log(1+Y4_3./Y4_2);
    Frac_YR_4=Frac_YR_4(2:jj);
    Intrg_G_4=[];
    Intrg_G_4(1)=(Y4_1(1)+Y4_1(2));
    Intrg_Y_4(1)=(Y4_3(1)+Y4_3(2));
    for i = 2:(jj-1)
        Intrg_G_4(i)=Intrg_G_4(i-1)+Y4_1(i)+Y4_1(i+1);
        Intrg_Y_4(i)=Intrg_Y_4(i-1)+Y4_3(i)+Y4_3(i+1);
    end
    Intrg_G_4=Intrg_G_4.*tstep./2;
    Intrg_Y_4=Intrg_Y_4.*tstep./2;
  
    %% plot frac vs integral
    Intrg_G=[Intrg_G_0;Intrg_G_1;Intrg_G_2;Intrg_G_3;Intrg_G_4];
    Intrg_Y=[Intrg_Y_0;Intrg_Y_1;Intrg_Y_2;Intrg_Y_3;Intrg_Y_4];
    Frac_YR=[Frac_YR_0';Frac_YR_1';Frac_YR_2';Frac_YR_3';Frac_YR_4'];
    Colo=['b','r','g','k','y'];
    figure;
%     subplot(1,3,1);
    for i = 1:5
       plot(Intrg_G(i,:)+Intrg_Y(i,:),Frac_YR(i,:),Colo(i),'LineWidth',2);
       hold on;
    end
    legend({'eta0','Ga Rr','Gr Ra','Ga Ra','Gr Rr'},'Fontsize',15);
    set(gca,'LineWidth',2,'Fontsize',20);
    xlabel('Inegral of Normalized G+Y Population','Fontsize',25);
    ylabel('[Y]/[R]','Fontsize',25);
    title(['Comparison between 5 Scenarios,dilution =',num2str(D)],'Fontsize',25);
       axis([0 40 0 2]);
%     axis([0 30 1e-10 4e9]);
%     
%     subplot(1,3,2);   
%     plot(Intrg_G_0+Intrg_Y_0,Frac_YR_0,'b',Intrg_G_1+Intrg_Y_1,Frac_YR_1,'r',Intrg_G_3+Intrg_Y_3,Frac_YR_3,'k','LineWidth',2.5);hold on;
%     % xlim([0 5E10]);
%     legend({'eta0','Ga Rr','Ga Ra'},'Fontsize',15);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     xlabel('Inegral of Normalized G+Y','Fontsize',25);
%     ylabel('[Y]/[R]','Fontsize',25);
%     title(['G Activation Conditions,[Cm]= ',num2str(A)],'Fontsize',25);
%     axis([0 35 0 1.0]);

% %     subplot(1,3,3);
%     plot(Intrg_G_0,Frac_YR_0,'b',Intrg_G_2,Frac_YR_2,'g',Intrg_G_4,Frac_YR_4,'y','LineWidth',2.5);hold on;
%     xlim([0 25]);
%     legend({'eta0','Gr Ra','Gr Rr'},'Fontsize',15);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     xlabel('Inegral of Normalized Donor Population','Fontsize',25);
%     ylabel('[Y]/[R]','Fontsize',25);
%     title(['G Repression Conditions,[Cm]= ',num2str(A)],'Fontsize',25);

    
      k(j)=getframe(gcf)
end
 movie2gif(k,'dilu_norm_total_test.gif','LoopCount',4,'DelayTime',1);

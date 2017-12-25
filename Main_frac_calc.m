%{ 
-------------main function----------
This is used in this set of function.
This is used for comparison among 5 systems.
Generate figures of ranking heatmap, bar and curve, and population
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
D=0.0;

%segregation error (don't change)
kappa=0.0;

% antibiotics type 
Atype = 'kan';
% A = 0; 
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


tspan= 0.1:0.1:60;
initConc = [1E-3*Nm; 1E-3*Nm; 0];
etaC_prime = 0.0375;
eta0 = (etaC_prime * mu_Y_max / Nm);
%% do integral
%{ 
Idea to generate parameter: to keep the eta_C integral as the same as
eta_0, basically can do double integral to the term, or seperately
calculate as followed\
%}

%This part should be commented if explore the effect of K1(K2)
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
% % KK=0.01:0.01:0.05;
% % kk=length(KK);
DD=0.0:0.005:0.2;
dd=length(DD);
AA=0.0:0.05:2.0;
aa=length(AA);

for i = 1:dd
      D=DD(i);
      for j = 1: aa
          A=AA(j);
           Version = 0;
           [GY_frac_0,Y_frac_0,mu_eff_update0,etaGR0,etaYR0]=Cell_frac(D,kappa,Version,[],[],[]);
           y0(j)=GY_frac_0;
           y00(j)=Y_frac_0;
            
            % G a R r
            Version = 1;
            [GY_frac_1,Y_frac_1,mu_eff_update1,etaGR1,etaYR1]=Cell_frac(D,kappa,Version,q_GA,q_RR,q_YA);
            y1(j)=GY_frac_1;
            y11(j)=Y_frac_1;
           
            %G r R a
            Version = 2;
            [GY_frac_2,Y_frac_2,mu_eff_update2,etaGR2,etaYR2]=Cell_frac(D,kappa,Version,q_GR,q_RA,q_YR);
            y2(j)=GY_frac_2;
            y22(j)=Y_frac_2;
                        
            %G a R a
            Version = 3;
            [GY_frac_3,Y_frac_3,mu_eff_update3,etaGR3,etaYR3]=Cell_frac(D,kappa,Version,q_GA,q_RA,q_YA);
            y3(j)=GY_frac_3;
            y33(j)=Y_frac_3;

            %G r R r
            Version = 4;
            [GY_frac_4,Y_frac_4,mu_eff_update4,etaGR4,etaYR4]=Cell_frac(D,kappa,Version,q_GR,q_RR,q_YR);
            y4(j)=GY_frac_4;
            y44(j)=Y_frac_4;
       
            % eta_total
            xx = etaYR1+etaGR1;
            yy =  etaYR2+etaGR2;
            zz = etaYR3+etaGR3;
            ww =  etaYR4+etaGR4;
      end
end

%% imagesc heatmap
%x is the variable in the second loop, z is the variable in the first loop
Popheat(AA,DD,'[Kan]','Dilution rate/Hr','GY',y0,y1,y2,y3,y4);
Popheat(AA,DD,'[Cm]','Dilution rate/Hr','GY',y00,y11,y22,y33,y44);

%% only change one parameter. Under different antibitics concentration
% DD=[0.0 0.05 0.1 0.15 0.2];
% dd=length(DD);
% AA=0.0:0.5:2.0;
% aa=length(AA);

% for j = 1: aa
%     if j ==1
%         Atype = 'none';
%         A=0.0;
%         Version = 0;
%         [GY_frac_0,Y_frac_0,mu_eff_update0,etaGR0,etaYR0]=Cell_frac(D,kappa,Version,[],[],[]);
%         y0(j)=GY_frac_0;
%         y00(j)=Y_frac_0;
%             
%             % G a R r
%             Version = 1;
%             [GY_frac_1,Y_frac_1,mu_eff_update1,etaGR1,etaYR1]=Cell_frac(D,kappa,Version,q_GA,q_RR,q_YA);
%             y1(j)=GY_frac_1;
%             y11(j)=Y_frac_1;
%            
%             %G r R a
%             Version = 2;
%             [GY_frac_2,Y_frac_2,mu_eff_update2,etaGR2,etaYR2]=Cell_frac(D,kappa,Version,q_GR,q_RA,q_YR);
%             y2(j)=GY_frac_2;
%             y22(j)=Y_frac_2;
%                         
%             %G a R a
%             Version = 3;
%             [GY_frac_3,Y_frac_3,mu_eff_update3,etaGR3,etaYR3]=Cell_frac(D,kappa,Version,q_GA,q_RA,q_YA);
%             y3(j)=GY_frac_3;
%             y33(j)=Y_frac_3;
% 
%             %G r R r
%             Version = 4;
%             [GY_frac_4,Y_frac_4,mu_eff_update4,etaGR4,etaYR4]=Cell_frac(D,kappa,Version,q_GR,q_RR,q_YR);
%             y4(j)=GY_frac_4;
%             y44(j)=Y_frac_4;
%        
%             % eta_total
%             xx = etaYR1+etaGR1;
%             yy =  etaYR2+etaGR2;
%             zz = etaYR3+etaGR3;
%             ww =  etaYR4+etaGR4;
%             
%     elseif j>1      
%           A=AA(j);
%           Atype = 'kan';
%             %constant etaC with dilution D
%             Version = 0;
%         [GY_frac_0,Y_frac_0,mu_eff_update0,etaGR0,etaYR0]=Cell_frac(D,kappa,Version,[],[],[]);
%         y0(j)=GY_frac_0;
%         y00(j)=Y_frac_0;
%             
%             % G a R r
%             Version = 1;
%             [GY_frac_1,Y_frac_1,mu_eff_update1,etaGR1,etaYR1]=Cell_frac(D,kappa,Version,q_GA,q_RR,q_YA);
%             y1(j)=GY_frac_1;
%             y11(j)=Y_frac_1;
%            
%             %G r R a
%             Version = 2;
%             [GY_frac_2,Y_frac_2,mu_eff_update2,etaGR2,etaYR2]=Cell_frac(D,kappa,Version,q_GR,q_RA,q_YR);
%             y2(j)=GY_frac_2;
%             y22(j)=Y_frac_2;
%                         
%             %G a R a
%             Version = 3;
%             [GY_frac_3,Y_frac_3,mu_eff_update3,etaGR3,etaYR3]=Cell_frac(D,kappa,Version,q_GA,q_RA,q_YA);
%             y3(j)=GY_frac_3;
%             y33(j)=Y_frac_3;
% 
%             %G r R r
%             Version = 4;
%             [GY_frac_4,Y_frac_4,mu_eff_update4,etaGR4,etaYR4]=Cell_frac(D,kappa,Version,q_GR,q_RR,q_YR);
%             y4(j)=GY_frac_4;
%             y44(j)=Y_frac_4;
%        
%             % eta_total
%             xx = etaYR1+etaGR1;
%             yy =  etaYR2+etaGR2;
%             zz = etaYR3+etaGR3;
%             ww =  etaYR4+etaGR4;
%     end
% end
%% under different dilution 
% Atype = 'none';
% A=0.0;
% DD= 0.0:0.01:0.15;
% dd=length(DD);
%  
% for j = 1: dd
%           D=DD(j);
%             %constant etaC with dilution D
%            Version = 0;
%            [GY_frac_0,Y_frac_0,mu_eff_update0,etaGR0,etaYR0]=Cell_frac(D,kappa,Version,[],[],[]);
%            y0(j)=GY_frac_0;
%            y00(j)=Y_frac_0;
%             
%             % G a R r
%             Version = 1;
%             [GY_frac_1,Y_frac_1,mu_eff_update1,etaGR1,etaYR1]=Cell_frac(D,kappa,Version,q_GA,q_RR,q_YA);
%             y1(j)=GY_frac_1;
%             y11(j)=Y_frac_1;
%            
%             %G r R a
%             Version = 2;
%             [GY_frac_2,Y_frac_2,mu_eff_update2,etaGR2,etaYR2]=Cell_frac(D,kappa,Version,q_GR,q_RA,q_YR);
%             y2(j)=GY_frac_2;
%             y22(j)=Y_frac_2;
%                         
%             %G a R a
%             Version = 3;
%             [GY_frac_3,Y_frac_3,mu_eff_update3,etaGR3,etaYR3]=Cell_frac(D,kappa,Version,q_GA,q_RA,q_YA);
%             y3(j)=GY_frac_3;
%             y33(j)=Y_frac_3;
% 
%             %G r R r
%             Version = 4;
%             [GY_frac_4,Y_frac_4,mu_eff_update4,etaGR4,etaYR4]=Cell_frac(D,kappa,Version,q_GR,q_RR,q_YR);
%             y4(j)=GY_frac_4;
%             y44(j)=Y_frac_4;
%        
%             % eta_total
%             xx = etaYR1+etaGR1;
%             yy =  etaYR2+etaGR2;
%             zz = etaYR3+etaGR3;
%             ww =  etaYR4+etaGR4;
% 
% end
%% plot dilution curves
% Poprank(DD,'Dilution rate/Hr','GY',y0,y1,y2,y3,y4);
% Poprank(DD,'Dilution rate/Hr','Y',y00,y11,y22,y33,y44);

%% plot ranking curves
% % Poprank(AA,'[Kan]','GY',y0,y1,y2,y3,y4);
% % Poprank(AA,'[Kan]','Y',y00,y11,y22,y33,y44);
%% plot bar
% % Popbar(AA,'[Kan]','GY',y0,y1,y2,y3,y4);
% % Popbar(AA,'[Kan]','Y',y00,y11,y22,y33,y44);

%% plot All figures
% % Popfig(tv,FvS,Fv0,Fv1,Fv2,Fv3,Fv4,xx,yy,zz,ww,etaGR1,etaGR2,etaGR3,etaGR4,etaYR1,etaYR2,etaYR3,etaYR4)

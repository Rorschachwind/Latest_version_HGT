function [GY_frac,Y_frac,mu_eff,etaGR,etaYR]=Cell_frac(D,kappa,Version,q_G,q_R,q_Y)
global tspan initConc
[tv,Fv]=ode45(@(t,Y) fun_GRY_Hill_D(t,Y,D,kappa,Version,q_G,q_R,q_Y),tspan,initConc);
Y_1=Fv(:,1);
Y_2=Fv(:,2);
Y_3=Fv(:,3);
ll=length(Y_1);
GY_frac=(Y_1(ll)+Y_3(ll))/(Y_1(ll)+Y_2(ll)+Y_3(ll));
Y_frac=Y_3(ll)/(Y_1(ll)+Y_2(ll)+Y_3(ll));

[mu_eff,etaGR,etaYR]=calcE(Fv,Version,q_G,q_R,q_Y);
return

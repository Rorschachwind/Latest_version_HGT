function dydt= fun_GRY_Hill_D(t,Y,D,kappa,version,q1,q2,q3)
%This function is the ODEs used in model

g=Y(1,1);
r=Y(2,1);
y=Y(3,1);

kappa=0;
[mu_eff,etaGR,etaYR]=fun_mu_Hill(g,r,y,version,q1,q2,q3);

mu_G_eff=mu_eff(1);
mu_R_eff=mu_eff(2);
mu_Y_eff=mu_eff(3);

dydt(1,1)=mu_G_eff*g - D*g - kappa*g;
dydt(2,1)=mu_R_eff*r - etaGR*r*g - etaYR*r*y -D*r + kappa*(g+y);
dydt(3,1)=mu_Y_eff*y + etaGR*r*g + etaYR*r*y - D*y -kappa*y;
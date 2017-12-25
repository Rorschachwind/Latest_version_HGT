function dydt = fun_GRY_Hill_D_update(t,Y,D,version,Fv)
%This function is used for update of the ODEs, where X is the multiply
%factor.
%This function is not used in the new version of modeling

X =calcX(Fv,version);

g=Y(1,1);
r=Y(2,1);
y=Y(3,1);

[mu_eff,etaGR,etaYR]=fun_mu_Hill_update(g,r,y,version,X);

mu_G_eff=mu_eff(1);
mu_R_eff=mu_eff(2);
mu_Y_eff=mu_eff(3);

dydt(1,1)=mu_G_eff*g - D*g;
dydt(2,1)=mu_R_eff*r - etaGR*r*g - etaYR*r*y -D*r;
dydt(3,1)=mu_Y_eff*y + etaGR*r*g + etaYR*r*y - D*y;



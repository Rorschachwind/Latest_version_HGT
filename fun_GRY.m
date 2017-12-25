function dydt= fun_GRY(t,Y,D,version)
%This function is the original version of the model
%Should pay attention to mu_R_max to change with antibiotics
global eta0 Nm

g=Y(1,1);
r=Y(2,1);
y=Y(3,1);

mu_G_max=0.33;
mu_Y_max=0.31;
mu_R_max=0.32;%kan=0.28

mu_G_eff=mu_G_max*(1-(g+r+y)/Nm);
mu_R_eff=mu_R_max*(1-(g+r+y)/Nm);
mu_Y_eff=mu_Y_max*(1-(g+r+y)/Nm);

dydt(1,1)=mu_G_eff*g - D*g;
dydt(2,1)=mu_R_eff*r - eta0*r*g - eta0*r*y -D*r;
dydt(3,1)=mu_Y_eff*y + eta0*r*g + eta0*r*y - D*y;
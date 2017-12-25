function X =calcX(Fv,version)
%This function is used to calculate the multiply factor in equations
%This function is not used in new version

global K1 K2 alpha1 alpha2 beta1 beta2 n m Nm mu_G_max mu_R_max mu_Y_max eta0 A Atype

% 1:G a R r
% 2: G r R a
% 3: G a R a
% 4: G r R r
% G activate R repress
mu_vec=fun_mu(A,Atype);
% mu_vec=[mu_G_max mu_R_max mu_Y_max];
v1=Fv(:,1);
v2=Fv(:,2);
v3=Fv(:,3);

mu_G_eff = mu_vec(1)*(1-(v1+v2+v3)/Nm);
mu_R_eff = mu_vec(2)*(1-(v1+v2+v3)/Nm);
mu_Y_eff = mu_vec(3)*(1-(v1+v2+v3)/Nm);

if version == 1
    Hill_G = alpha1 + alpha2 * mu_G_eff.^n./(K1^n+mu_G_eff.^n);
    Hill_R = beta1 + beta2 * K2^m./(K2^m+mu_R_eff.^m);
    Hill_Y = alpha1 + alpha2 * mu_Y_eff.^n./(K1^n+mu_Y_eff.^n);

elseif version == 2
    Hill_G = alpha1 + alpha2 * K1^n./(K1^n+mu_G_eff.^n);
    Hill_R = beta1 + beta2 * mu_R_eff.^m./(K2^m+mu_R_eff.^m);
    Hill_Y = alpha1 + alpha2 * mu_Y_eff.^n./(K1^n+mu_Y_eff.^n);

elseif version == 3
    Hill_G = alpha1 + alpha2 * mu_G_eff.^n./(K1^n+mu_G_eff.^n);
    Hill_R = beta1 + beta2 * mu_R_eff.^m./(K2^m+mu_R_eff.^m);
    Hill_Y = alpha1 + alpha2 * mu_Y_eff.^n./(K1^n+mu_Y_eff.^n);
    
elseif version == 4
    Hill_G = alpha1 + alpha2 * K1^n./(K1^n+mu_G_eff.^n);
    Hill_R = beta1 + beta2 * K2^m./(K2^m+mu_R_eff.^m);
    Hill_Y = alpha1 + alpha2 * mu_Y_eff.^n./(K1^n+mu_Y_eff.^n);
end

eta_Hill_GR = eta0.*Hill_G.*Hill_R;
eta_Hill_YR = eta0.*Hill_Y.*Hill_R;

X=sqrt(eta0/mean(eta_Hill_GR));
% A1=alpha1*X;
% A2=alpha2*X;
% B1=beta1*X;
% B2=beta2*X;
return




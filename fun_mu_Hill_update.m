function [mu_eff,etaGR,etaYR]=fun_mu_Hill_update(g,r,y,version,X)

%This function is used for the update of parmeters
%This function is not used for new version modeling
% 1:G a R r
% 2: G r R a
% 3: G a R a
% 4: G r R r
%G activte R repress
global K1 K2 alpha1 alpha2 beta1 beta2 n m Nm mu_G_max mu_R_max mu_Y_max eta0 A Atype

mu_vec=fun_mu(A,Atype);
mu_eff = [mu_vec(1)*(1-(g+r+y)/Nm) mu_vec(2)*(1-(g+r+y)/Nm) mu_vec(3)*(1-(g+r+y)/Nm)];

mu_G_eff=mu_eff(1);
mu_R_eff=mu_eff(2);
mu_Y_eff=mu_eff(3);

if version == 1
    Hill_G = alpha1 * X + alpha2 * X * mu_G_eff^n/(K1^n+mu_G_eff^n);
    Hill_R = beta1 * X + beta2 * X * K2^m/(K2^m+mu_R_eff^m);
    Hill_Y = alpha1 * X + alpha2 * X * mu_Y_eff^n/(K1^n+mu_Y_eff^n);
elseif version == 2
    Hill_G = alpha1 * X + alpha2 * X * K1^n/(K1^n+mu_G_eff^n);
    Hill_R = beta1 * X + beta2 * X * mu_R_eff.^m/(K2^m+mu_R_eff^m);
    Hill_Y = alpha1 * X + alpha2 * X * mu_Y_eff^n/(K1^n+mu_Y_eff^n);
elseif version == 3
    Hill_G = alpha1 * X + alpha2 * X * mu_G_eff^n/(K1^n+mu_G_eff^n);
    Hill_R = beta1 * X + beta2 * X * mu_R_eff.^m/(K2^m+mu_R_eff^m);
    Hill_Y = alpha1 * X + alpha2 * X * mu_Y_eff^n/(K1^n+mu_Y_eff^n);
elseif version == 4
    Hill_G = alpha1 * X + alpha2 * X * K1^n/(K1^n+mu_G_eff^n);
    Hill_R = beta1 * X + beta2 * X * K2^m/(K2^m+mu_R_eff^m);
    Hill_Y = alpha1 * X + alpha2 * X * mu_Y_eff^n/(K1^n+mu_Y_eff^n);
end


etaGR = eta0*Hill_G*Hill_R;
etaYR = eta0*Hill_Y*Hill_R;

return
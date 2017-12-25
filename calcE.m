function [mu_eff_update,etaGR,etaYR]= calcE (Fv,version,q1,q2,q3)
%This function is used to calculate and display the resultds for etaGR and
%etaYR,also mu_eff

global K1 K2 alpha1 alpha2 beta1 beta2 n m Nm mu_G_max mu_R_max mu_Y_max eta0 A Atype

% G activate R repress
% mu_vec=[mu_G_max mu_R_max mu_Y_max];
mu_vec=fun_mu(A,Atype);
v1=Fv(:,1);
v2=Fv(:,2);
v3=Fv(:,3);

mu_G_eff = mu_vec(1)*(1-(v1+v2+v3)/Nm);
mu_R_eff = mu_vec(2)*(1-(v1+v2+v3)/Nm);
mu_Y_eff = mu_vec(3)*(1-(v1+v2+v3)/Nm);


if version == 1
    Hill_G = alpha1 /q1 + alpha2 /q1 * mu_G_eff.^n./(K1^n+mu_G_eff.^n);
    Hill_R = beta1 /q2 + beta2 /q2 * K2^m./(K2^m+mu_R_eff.^m);
    Hill_Y = alpha1/q3  + alpha2/q3  * mu_Y_eff.^n./(K1^n+mu_Y_eff.^n);
elseif version == 2
    Hill_G = alpha1 /q1 + alpha2 /q1 * K1^n./(K1^n+mu_G_eff.^n);
    Hill_R = beta1 /q2 + beta2 /q2 * mu_R_eff.^m./(K2^m+mu_R_eff.^m);
    Hill_Y = alpha1/q3  + alpha2/q3  * K1.^n./(K1^n+mu_Y_eff.^n);
elseif version == 3
    Hill_G = alpha1 /q1 + alpha2 /q1 * mu_G_eff.^n./(K1^n+mu_G_eff.^n);
    Hill_R = beta1 /q2 + beta2 /q2 * mu_R_eff.^m./(K2^m+mu_R_eff.^m);
    Hill_Y = alpha1/q3  + alpha2/q3  * mu_Y_eff.^n./(K1^n+mu_Y_eff.^n);
elseif version == 4
    Hill_G = alpha1 /q1 + alpha2 /q1 * K1^n./(K1^n+mu_G_eff.^n);
    Hill_R = beta1 /q2 + beta2 /q2 * K2^m./(K2^m+mu_R_eff.^m);
    Hill_Y = alpha1/q3  + alpha2/q3  * K1.^n./(K1^n+mu_Y_eff.^n);
elseif version == 0
    Hill_G=1;
    Hill_R=1;
    Hill_Y=1;
end

mu_eff_update=[mu_G_eff mu_R_eff mu_Y_eff];
etaGR = eta0 .* Hill_G .* Hill_R;
etaYR = eta0 .* Hill_Y .* Hill_R;

return
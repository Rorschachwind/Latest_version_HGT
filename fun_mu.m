function mus = fun_mu(A, Atype)

%This funtion is used to generate mu_max part

if isequal(lower(Atype), 'cm')
    % G is sensitive to Cm
    % R is resistant to Cm
    % Y is resistant to Cm
    
    mu_G_max = 0.33;
    K_G = 1.92;
    n_G = 2.19;
    
%     mu_R_max_Cm = 0.28;
     mu_R_max_Cm = 0.32;
    m_R_Cm = 2.17E-4;
    
    mu_Y_max = 0.31;
    m_Y_c = 2.31E-4;
    
    mus(1) = mu_G_max * K_G^n_G / (K_G^n_G + A^n_G); %G; cm % when A = 0.5, = 0.3135
    mus(2) = mu_R_max_Cm - m_R_Cm * A; %R; cm % when A = 0.5, =0.3199
    mus(3) = mu_Y_max - m_Y_c * A; %Y; cm % when A = 0.5, =0.3099
    
elseif isequal(lower(Atype), 'kan')
    % G is resistant to Kan
    % R is sensitive to Kan
    % Y is resistant to Kan
    
    mu_G_max = 0.33;
    m_G = 1.86E-5;
    
    mu_R_max_Kan = 0.28;
    K_R = 2.01;
    n_R = 6.45;
    
    mu_Y_max = 0.31;
    m_Y_k = 2.66E-5;
    
    mus(1) = mu_G_max - m_G * A;
    mus(2) = mu_R_max_Kan * K_R^n_R / (K_R^n_R + A^n_R);
    mus(3) = mu_Y_max + m_Y_k * A;
    
elseif isequal(lower(Atype), 'none')
    
    mu_G_max = 0.33;
    m_G = 1.86E-5;
    
     mu_R_max_Cm = 0.32;
%     mu_R_max_Cm = 0.28;
    m_R_Cm = 2.17E-4;
    
    mu_Y_max = 0.31;
    m_Y_k = 2.66E-5;
    
    mus(1) = mu_G_max - m_G * A; %G; kan
    mus(2) = mu_R_max_Cm - m_R_Cm * A; %R; cm
    mus(3) = mu_Y_max + m_Y_k * A; %Y; kan
    
elseif isequal(lower(Atype), 'both')
    
    mu_G_max = 0.33;
    K_G = 1.92;
    n_G = 2.19;
    
    mu_R_max_Kan = 0.28;
    K_R = 2.01;
    n_R = 6.45;
    
    mu_Y_max = 0.31;
    m_Y_k = 2.66E-5;
    
    mus(1) = mu_G_max * K_G^n_G / (K_G^n_G + A^n_G); %G; cm
    mus(2) = mu_R_max_Kan * K_R^n_R / (K_R^n_R + A^n_R); %R; kan
    mus(3) = mu_Y_max + m_Y_k * A; %Y; kan
    
else
    error('Error: Atype should be either "None", "Cm", "Kan", or "Both"')
    
end

%% G
% mu_G_max = 0.33;
% K_G = 1.92;
% n_G = 2.19;
% m_G = 1.86E-5;
% m_G_Cm = 6.57E-6;
% mu_G_max_c = 0.15;

% mu_G_Cm = mu_G_max * K_G^n_G / (K_G^n_G + A^n_G);
% mu_G_Kan = mu_G_max - m_G * A;
% mu_G_KanCm = mu_G_max_c - m_G_Cm;

%% R
% mu_R_max_Cm = 0.32;
% mu_R_Cm = 2.17E-4;
% mu_R_max_Kan = 0.28;
% K_R = 2.01;
% n_R = 6.45;
% mu_R_Cm = mu_R_max_Cm - m_R_Cm * A;
% mu_R_Kan = mu_R_max_Kan * K_R^n_R / (K_R^n_R + A^n_R);

%% Y
% mu_Y_max = 0.31;
% m_Y_c = 2.31E-4;
% m_Y_k = 2.66E-5;
% mu_Y_Cm = mu_Y_max - m_Y_c * A;
% mu_Y_Kan = mu_Y_max + m_Y_k * A;



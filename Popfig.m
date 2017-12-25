function f= Popfig(tv,FvS,Fv0,Fv1,Fv2,Fv3,Fv4,xx,yy,zz,ww,eta1,eta2,eta3,eta4,eta11,eta22,eta33,eta44)
%{
This function is used for generate figures in this problem.
%}

global eta0


    % G activate R repress eta_C constant no dilution (Original)
    fracS_G=FvS(:,1)./(FvS(:,1)+FvS(:,2)+FvS(:,3));
    fracS_R=FvS(:,2)./(FvS(:,1)+FvS(:,2)+FvS(:,3));
    fracS_Y=FvS(:,3)./(FvS(:,1)+FvS(:,2)+FvS(:,3));
    % G activate R repress eta_C constant
    frac0_G=Fv0(:,1)./(Fv0(:,1)+Fv0(:,2)+Fv0(:,3));
    frac0_R=Fv0(:,2)./(Fv0(:,1)+Fv0(:,2)+Fv0(:,3));
    frac0_Y=Fv0(:,3)./(Fv0(:,1)+Fv0(:,2)+Fv0(:,3));
    % G activate R repress eta_C function
    frac1_G=Fv1(:,1)./(Fv1(:,1)+Fv1(:,2)+Fv1(:,3));
    frac1_R=Fv1(:,2)./(Fv1(:,1)+Fv1(:,2)+Fv1(:,3));
    frac1_Y=Fv1(:,3)./(Fv1(:,1)+Fv1(:,2)+Fv1(:,3));
    %G repress R activate eta_C function
    frac2_G=Fv2(:,1)./(Fv2(:,1)+Fv2(:,2)+Fv2(:,3));
    frac2_R=Fv2(:,2)./(Fv2(:,1)+Fv2(:,2)+Fv2(:,3));
    frac2_Y=Fv2(:,3)./(Fv2(:,1)+Fv2(:,2)+Fv2(:,3));
    %G activate R activate eta_C function
    frac3_G=Fv3(:,1)./(Fv3(:,1)+Fv3(:,2)+Fv3(:,3));
    frac3_R=Fv3(:,2)./(Fv3(:,1)+Fv3(:,2)+Fv3(:,3));
    frac3_Y=Fv3(:,3)./(Fv3(:,1)+Fv3(:,2)+Fv3(:,3));
    %G repress R repress eta_C function
    frac4_G=Fv4(:,1)./(Fv4(:,1)+Fv4(:,2)+Fv4(:,3));
    frac4_R=Fv4(:,2)./(Fv4(:,1)+Fv4(:,2)+Fv4(:,3));
    frac4_Y=Fv4(:,3)./(Fv4(:,1)+Fv4(:,2)+Fv4(:,3));

    % figure1 3 popu frac
    figure;
    subplot(1,5,1);
    plot(tv,frac0_G,'g',tv,frac0_R,'r',tv,frac0_Y,'y',tv,fracS_G,'g:',tv,fracS_R,'r:',tv,fracS_Y,'y:','LineWidth',2);
    xlabel('Time/min','Fontsize',20);
    ylabel('Fraction of cells','Fontsize',20);
    ylim([0 1]);
    set(gca,'LineWidth',2,'Fontsize',20);
    title('\eta_C = \eta_0,G: activation|R: repression','Fontsize',20);
    legend('G','R','Y');
    subplot(1,5,2);
    plot(tv,frac1_G,'g',tv,frac1_R,'r',tv,frac1_Y,'y','LineWidth',2);
    hold on;
    plot(tv,frac0_G,'--',tv,frac0_R,'--',tv,frac0_Y,'--','LineWidth',1.5);
    ylim([0 1]);
    xlabel('Time/min','Fontsize',20);
    ylabel('Fraction of cells','Fontsize',20);
    set(gca,'LineWidth',2,'Fontsize',20);
    title('\eta_C = \eta(t),G: activation|R: repression','Fontsize',20);
    % legend('G','R','Y');
    subplot(1,5,3);
    plot(tv,frac2_G,'g',tv,frac2_R,'r',tv,frac2_Y,'y','LineWidth',2);
    hold on;
    plot(tv,frac0_G,'--',tv,frac0_R,'--',tv,frac0_Y,'--','LineWidth',1.5);
    ylim([0 1]);
    xlabel('Time/min','Fontsize',20);
    ylabel('Fraction of cells','Fontsize',20);
    set(gca,'LineWidth',2,'Fontsize',20);
    title('\eta_C = \eta(t),G: repression|R: activation','Fontsize',20);
    % legend('G','R','Y');
    subplot(1,5,4);
    plot(tv,frac3_G,'g',tv,frac3_R,'r',tv,frac3_Y,'y','LineWidth',2);
    hold on;
    plot(tv,frac0_G,'--',tv,frac0_R,'--',tv,frac0_Y,'--','LineWidth',1.5);
    ylim([0 1]);
    xlabel('Time/min','Fontsize',20);
    ylabel('Fraction of cells','Fontsize',20);
    set(gca,'LineWidth',2,'Fontsize',20);
    title('\eta_C = \eta(t),G: activation|R: activation','Fontsize',20);
    % legend('G','R','Y');
    subplot(1,5,5);
    plot(tv,frac4_G,'g',tv,frac4_R,'r',tv,frac4_Y,'y','LineWidth',2);
    hold on;
    plot(tv,frac0_G,'--',tv,frac0_R,'--',tv,frac0_Y,'--','LineWidth',1.5);
    ylim([0 1]);
    xlabel('Time/min','Fontsize',20);
    ylabel('Fraction of cells','Fontsize',20);
    set(gca,'LineWidth',2,'Fontsize',20);
    title('\eta_C = \eta(t),G: repression|R: repression','Fontsize',20);
    suptitle('Fraction of 3 Population');
    h=suptitle('Fraction of 3 Population');
    set(h,'Fontsize',25);

    % figure2 fraction of with/without plasmid
    figure;
    subplot(1,5,1);
    plot(tv,frac0_G+frac0_Y,'b',tv,frac0_R,'r',tv,fracS_G+fracS_Y,'b:',tv,fracS_R,'r:','LineWidth',2);
    ylim([0 1]);
    xlabel('Time/min','Fontsize',20);
    ylabel('Fraction of cells','Fontsize',20);
    set(gca,'LineWidth',2,'Fontsize',20);
    title('\eta_C = \eta_0,G: activation|R: repression','Fontsize',20);
    legend('G+Y','R');
    subplot(1,5,2);
    plot(tv,frac1_G+frac1_Y,'b',tv,frac1_R,'r','LineWidth',2);
    hold on;
    plot(tv,frac0_G+frac0_Y,'--',tv,frac0_R,'--','LineWidth',1.5);
    ylim([0 1]);
    xlabel('Time/min','Fontsize',20);
    ylabel('Fraction of cells','Fontsize',20);
    set(gca,'LineWidth',2,'Fontsize',20);
    title('\eta_C = \eta(t),G: activation|R: repression','Fontsize',20);
    % legend('G+Y','R');
    subplot(1,5,3);
    plot(tv,frac2_G+frac2_Y,'b',tv,frac2_R,'r','LineWidth',2);
    hold on;
    plot(tv,frac0_G+frac0_Y,'--',tv,frac0_R,'--','LineWidth',1.5);
    ylim([0 1]);
    xlabel('Time/min','Fontsize',20);
    ylabel('Fraction of cells','Fontsize',20);
    set(gca,'LineWidth',2,'Fontsize',20);
    title('\eta_C = \eta(t),G: repression|R: activation','Fontsize',20);
    % legend('G+Y','R');
    subplot(1,5,4);
    plot(tv,frac3_G+frac3_Y,'b',tv,frac3_R,'r','LineWidth',2);
    hold on;
    plot(tv,frac0_G+frac0_Y,'--',tv,frac0_R,'--','LineWidth',1.5);
    ylim([0 1]);
    xlabel('Time/min','Fontsize',20);
    ylabel('Fraction of cells','Fontsize',20);
    set(gca,'LineWidth',2,'Fontsize',20);
    title('\eta_C = \eta(t),G: activation|R: activation','Fontsize',20);
    % legend('G+Y','R');
    subplot(1,5,5);
    plot(tv,frac4_G+frac4_Y,'b',tv,frac4_R,'r','LineWidth',2);
    hold on;
    plot(tv,frac0_G+frac0_Y,'--',tv,frac0_R,'--','LineWidth',1.5);
    ylim([0 1]);
    xlabel('Time/min','Fontsize',20);
    ylabel('Fraction of cells','Fontsize',20);
    set(gca,'LineWidth',2,'Fontsize',20);
    title('\eta_C = \eta(t),G: repression|R: repression','Fontsize',20);
    suptitle('Fraction of Plasmid-containing Population');
    h=suptitle('Fraction of Plasmid-containing Population');
    set(h,'Fontsize',25);

    % figure3 plot those 2 on the same figure
    figure;
    subplot(3,5,1);
    plot(tv,frac0_G,'g',tv,frac0_R,'r',tv,frac0_Y,'y',tv,fracS_G,'g:',tv,fracS_R,'r:',tv,fracS_Y,'y:','LineWidth',2);
    ylim([0 1]);
%     xlabel('Time/min','Fontsize',15);
    ylabel('Fraction of cells','Fontsize',15);
    set(gca,'LineWidth',2,'Fontsize',15);
    title('\eta_C = \eta_0,G: activation|R: repression','Fontsize',15);
    legend('G','R','Y');
    subplot(3,5,2);
    plot(tv,frac1_G,'g',tv,frac1_R,'r',tv,frac1_Y,'y','LineWidth',1.5);
    hold on;
    plot(tv,frac0_G,'--',tv,frac0_R,'--',tv,frac0_Y,'--','LineWidth',1.5);
    ylim([0 1]);
%     xlabel('Time/min','Fontsize',15);
%      ylabel('Fraction of cells','Fontsize',15);
    set(gca,'LineWidth',2,'Fontsize',15);
    title('\eta_C = \eta(t),G: activation|R: repression','Fontsize',15);
    % legend('G','R','Y');
    subplot(3,5,3);
    plot(tv,frac2_G,'g',tv,frac2_R,'r',tv,frac2_Y,'y','LineWidth',2);
    hold on;
    plot(tv,frac0_G,'--',tv,frac0_R,'--',tv,frac0_Y,'--','LineWidth',1.5);
    ylim([0 1]);
%     xlabel('Time/min','Fontsize',15);
%     ylabel('Fraction of cells','Fontsize',15);
    set(gca,'LineWidth',2,'Fontsize',15);
    title('\eta_C = \eta(t),G: repression|R: activation','Fontsize',15);
    % legend('G','R','Y');
    subplot(3,5,4);
    plot(tv,frac3_G,'g',tv,frac3_R,'r',tv,frac3_Y,'y','LineWidth',2);
    hold on;
    plot(tv,frac0_G,'--',tv,frac0_R,'--',tv,frac0_Y,'--','LineWidth',1.5);
    ylim([0 1]);
%     xlabel('Time/min','Fontsize',15);
%     ylabel('Fraction of cells','Fontsize',15);
    set(gca,'LineWidth',2,'Fontsize',15);
    title('\eta_C = \eta(t),G: activation|R: activation','Fontsize',15);
    % legend('G','R','Y');
    subplot(3,5,5);
    plot(tv,frac4_G,'g',tv,frac4_R,'r',tv,frac4_Y,'y','LineWidth',2);
    hold on;
    plot(tv,frac0_G,'--',tv,frac0_R,'--',tv,frac0_Y,'--','LineWidth',1.5);
    ylim([0 1]);
%     xlabel('Time/min','Fontsize',15);
%     ylabel('Fraction of cells','Fontsize',15);
    set(gca,'LineWidth',2,'Fontsize',15);
    title('\eta_C = \eta(t),G: repression|R: repression','Fontsize',15);
    % legend('G','R','Y');

    subplot(3,5,6);
    plot(tv,frac0_G+frac0_Y,'b',tv,frac0_R,'r',tv,fracS_G+fracS_Y,'b:',tv,fracS_R,'r:','LineWidth',2);
    ylim([0 1]);
%     xlabel('Time/min','Fontsize',15);
    ylabel('Fraction of cells','Fontsize',15);
    set(gca,'LineWidth',2,'Fontsize',15);
%     title('\eta_C = \eta_0,G: activation|R: repression','Fontsize',15);
    legend('G+Y','R');
    subplot(3,5,7);
    plot(tv,frac1_G+frac1_Y,'b',tv,frac1_R,'r','LineWidth',2);
    hold on;
    plot(tv,frac0_G+frac0_Y,'--',tv,frac0_R,'--','LineWidth',1.5);
    ylim([0 1]);
%     xlabel('Time/min','Fontsize',15);
%     ylabel('Fraction of cells','Fontsize',15);
    set(gca,'LineWidth',2,'Fontsize',15);
%     title('\eta_C = \eta(t),G: activation|R: repression','Fontsize',15);
    % legend('G+Y','R');
    subplot(3,5,8);
    plot(tv,frac2_G+frac2_Y,'b',tv,frac2_R,'r','LineWidth',2);
    hold on;
    plot(tv,frac0_G+frac0_Y,'--',tv,frac0_R,'--','LineWidth',1.5);
    ylim([0 1]);
%     xlabel('Time/min','Fontsize',15);
%     ylabel('Fraction of cells','Fontsize',15);
    set(gca,'LineWidth',2,'Fontsize',15);
%     title('\eta_C = \eta(t),G: repression|R: activation','Fontsize',15);
    % legend('G+Y','R');
    subplot(3,5,9);
    plot(tv,frac3_G+frac3_Y,'b',tv,frac3_R,'r','LineWidth',2);
    hold on;
    plot(tv,frac0_G+frac0_Y,'--',tv,frac0_R,'--','LineWidth',1.5);
    ylim([0 1]);
%     xlabel('Time/min','Fontsize',15);
%     ylabel('Fraction of cells','Fontsize',15);
    set(gca,'LineWidth',2,'Fontsize',15);
%     title('\eta_C = \eta(t),G: activation|R: activation','Fontsize',15);
    % legend('G+Y','R');
    subplot(3,5,10);
    plot(tv,frac4_G+frac4_Y,'b',tv,frac4_R,'r','LineWidth',2);
    hold on;
    plot(tv,frac0_G+frac0_Y,'--',tv,frac0_R,'--','LineWidth',1.5);
    ylim([0 1]);
%     xlabel('Time/min','Fontsize',);
%     ylabel('Fraction of cells','Fontsize',15);
    set(gca,'LineWidth',2,'Fontsize',15);
%     title('\eta_C = \eta(t),G: repression|R: repression','Fontsize',15);
%     suptitle('Fraction of Cells');
%     h=suptitle('Fraction of Cells');
%     set(h,'Fontsize',25);
%     
    l=length(Fv0(:,1));
    subplot(3,5,11);
    semilogy(tv,eta0*ones(1,l),'LineWidth',2);
    xlabel('Time/min','Fontsize',15);
    ylabel('\eta_0','Fontsize',15);
    set(gca,'LineWidth',2,'Fontsize',15);
     ylim([1e-15 1e-5]);
%     title('G: activation|R: repression','Fontsize',20);
%     ylim([0 1e-10]);
    legend('\eta_0');
    subplot(3,5,12)
    semilogy(tv, xx,tv,eta0*ones(1,l),'--','LineWidth',2);
    xlabel('Time/min','Fontsize',15);
%     ylabel('\eta_{total}','Fontsize',20);
    set(gca,'LineWidth',2,'Fontsize',15);
     ylim([1e-15 1e-5]);
%     title('G: activation|R: repression','Fontsize',20);
%     ylim([0 1e-10]);
%     legend('\eta_{total}','\eta_0');
    subplot(3,5,13)
    semilogy(tv, yy,tv,eta0*ones(1,l),'--','LineWidth',2)
    xlabel('Time/min','Fontsize',15);
%     ylabel('\eta_{total}','Fontsize',20);
%     title('G: repression|R: activation','Fontsize',20);
    set(gca,'LineWidth',2,'Fontsize',15);
     ylim([1e-15 1e-5]);
%     ylim([0,1e-10]);
%     legend('\eta_{total}','\eta_0');
    subplot(3,5,14)
    semilogy(tv, zz,tv,eta0*ones(1,l),'--','LineWidth',2)
    xlabel('Time/min','Fontsize',15);
%     ylabel('\eta_{total}','Fontsize',20);
%     title('G: activation|R: activation','Fontsize',20);
    set(gca,'LineWidth',2,'Fontsize',15);
     ylim([1e-15 1e-5]);
%     ylim([0 1e-10]);
%     legend('\eta_{total}','\eta_0');
    subplot(3,5,15)
    semilogy(tv, ww,tv,eta0*ones(1,l),'--','LineWidth',2)
    xlabel('Time/min','Fontsize',15);
%     ylabel('\eta_{total}','Fontsize',20);
%     title('G: repression|R: repression','Fontsize',20);
    set(gca,'LineWidth',2,'Fontsize',15);
     ylim([1e-15 1e-5]);
%     legend('\eta_{total}','\eta_0');
%     suptitle('\eta_{total} v. \eta_0');
%     h=suptitle('\eta_{total} v. \eta_0');
%     set(h,'Fontsize',25);
    
    
    
    
    
    
%     % figure4 eta_total & eta0
%     l=length(Fv0(:,1));
%     figure;
%     subplot(1,5,1);
%     plot(tv,eta0*ones(1,l),'LineWidth',2);
%     xlabel('Time/min','Fontsize',20);
%     ylabel('\eta_0','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
% %     title('G: activation|R: repression','Fontsize',20);
%     ylim([0 1e-10]);
%     legend('\eta_0');
%     subplot(1,5,2)
%     plot(tv, xx,tv,eta0*ones(1,l),'--','LineWidth',2);
%     xlabel('Time/min','Fontsize',20);
% %     ylabel('\eta_{total}','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
% %     title('G: activation|R: repression','Fontsize',20);
%     ylim([0 1e-10]);
% %     legend('\eta_{total}','\eta_0');
%     subplot(1,5,3)
%     plot(tv, yy,tv,eta0*ones(1,l),'--','LineWidth',2)
%     xlabel('Time/min','Fontsize',20);
% %     ylabel('\eta_{total}','Fontsize',20);
% %     title('G: repression|R: activation','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     ylim([0,1e-10]);
% %     legend('\eta_{total}','\eta_0');
%     subplot(1,5,4)
%     plot(tv, zz,tv,eta0*ones(1,l),'--','LineWidth',2)
%     xlabel('Time/min','Fontsize',20);
% %     ylabel('\eta_{total}','Fontsize',20);
% %     title('G: activation|R: activation','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     ylim([0 1e-10]);
% %     legend('\eta_{total}','\eta_0');
%     subplot(1,5,5)
%     plot(tv, ww,tv,eta0*ones(1,l),'--','LineWidth',2)
%     xlabel('Time/min','Fontsize',20);
% %     ylabel('\eta_{total}','Fontsize',20);
% %     title('G: repression|R: repression','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     ylim([0 1e-10]);
% %     legend('\eta_{total}','\eta_0');
%     suptitle('\eta_{total} v. \eta_0');
%     h=suptitle('\eta_{total} v. \eta_0');
%     set(h,'Fontsize',25);
% 
% 
% %% plot etaYR etaGR
% 
%     l=length(tv);
%     figure;
%     subplot(1,5,1);
%     plot(tv,eta0*ones(1,l),'LineWidth',2);
%     xlabel('Time/min','Fontsize',20);
%     ylabel('\eta_0','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     title('G: activation|R: repression','Fontsize',20);
%     ylim([0 1e-10]);
%     legend('\eta_0');
%     subplot(1,5,2)
%     plot(tv, eta1,tv,eta0*ones(1,l),'--','LineWidth',2);
%     xlabel('Time/min','Fontsize',20);
%     ylabel('\eta_{GR}','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     title('G: activation|R: repression','Fontsize',20);
%     ylim([0 1e-10]);
%     legend('\eta_{GR}','\eta_0');
%     subplot(1,5,3)
%     plot(tv, eta2,tv,eta0*ones(1,l),'--','LineWidth',2)
%     xlabel('Time/min','Fontsize',20);
%     ylabel('\eta_{GR}','Fontsize',20);
%     title('G: repression|R: activation','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     ylim([0,1e-10]);
%     legend('\eta_{GR}','\eta_0');
%     subplot(1,5,4)
%     plot(tv, eta3,tv,eta0*ones(1,l),'--','LineWidth',2)
%     xlabel('Time/min','Fontsize',20);
%     ylabel('\eta_{GR}','Fontsize',20);
%     title('G: activation|R: activation','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     ylim([0 1e-10]);
%     legend('\eta_{GR}','\eta_0');
%     subplot(1,5,5)
%     plot(tv, eta4,tv,eta0*ones(1,l),'--','LineWidth',2)
%     xlabel('Time/min','Fontsize',20);
%     ylabel('\eta_{GR}','Fontsize',20);
%     title('G: repression|R: repression','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     ylim([0 1e-10]);
%     legend('\eta_{GR}','\eta_0');
%     suptitle('\eta_{GR} v. \eta_0');
%     h=suptitle('\eta_{GR} v. \eta_0');
%     set(h,'Fontsize',25);
% 
%  %figure for eta-Y
%  
%  l=length(tv);
%     figure;
%     subplot(1,5,1);
%     plot(tv,eta0*ones(1,l),'LineWidth',2);
%     xlabel('Time/min','Fontsize',20);
%     ylabel('\eta_0','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     title('G: activation|R: repression','Fontsize',20);
%     ylim([0 1e-10]);
%     legend('\eta_0');
%     subplot(1,5,2)
%     plot(tv, eta11,tv,eta0*ones(1,l),'--','LineWidth',2);
%     xlabel('Time/min','Fontsize',20);
%     ylabel('\eta_{YR}','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     title('G: activation|R: repression','Fontsize',20);
%     ylim([0 1e-10]);
%     legend('\eta_{YR}','\eta_0');
%     subplot(1,5,3)
%     plot(tv, eta22,tv,eta0*ones(1,l),'--','LineWidth',2)
%     xlabel('Time/min','Fontsize',20);
%     ylabel('\eta_{YR}','Fontsize',20);
%     title('G: repression|R: activation','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     ylim([0,1e-10]);
%     legend('\eta_{YR}','\eta_0');
%     subplot(1,5,4)
%     plot(tv, eta33,tv,eta0*ones(1,l),'--','LineWidth',2)
%     xlabel('Time/min','Fontsize',20);
%     ylabel('\eta_{YR}','Fontsize',20);
%     title('G: activation|R: activation','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     ylim([0 1e-10]);
%     legend('\eta_{YR}','\eta_0');
%     subplot(1,5,5)
%     plot(tv, eta44,tv,eta0*ones(1,l),'--','LineWidth',2)
%     xlabel('Time/min','Fontsize',20);
%     ylabel('\eta_{YR}','Fontsize',20);
%     title('G: repression|R: repression','Fontsize',20);
%     set(gca,'LineWidth',2,'Fontsize',20);
%     ylim([0 1e-10]);
%     legend('\eta_{YR}','\eta_0');
%     suptitle('\eta_{YR} v. \eta_0');
%     h=suptitle('\eta_{YR} v. \eta_0');
%     set(h,'Fontsize',25);
% 
%  
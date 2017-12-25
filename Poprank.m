function f = Poprank(x,name,frac_type,y0,y1,y2,y3,y4)

% name can be: [Kan],[Cm],Dilution rate/Hr
YYY=[y0;y1;y2;y3;y4];
Name=["G: activation | R: repression","G: repression | R: activation","G: activation | R: activation","G: repression | R: repression"];

figure;
% seperate subplot
for i = 1 :4
    subplot(1,4,i);
    plot(x,YYY(i+1,:),x,YYY(1,:),'--','LineWidth',2);
    ylim([0.0,1.0]);
    set(gca,'LineWidth',2,'Fontsize',20);
    xlabel(name,'Fontsize',20);
    title(Name(i),'Fontsize',25);
    xlim([0 max(x)]);
end
subplot(1,4,1);
H=legend('\eta(t)','\eta_0');
h=suptitle(['Fraction of cells under different ',name]);
set(h,'Fontsize',25);

if frac_type == 'GY'
    ylabel('Fraction of G+Y cells','Fontsize',20);
elseif frac_type == 'Y'
    ylabel('Fraction of Y cells','Fontsize',20);
end

%plot on the same figure
figure;
for i = 1 :4
    plot(x,YYY(i+1,:),'LineWidth',2);
    ylim([0.0,1.0]);
    set(gca,'LineWidth',2,'Fontsize',20);
    xlim([0 max(x)]);
    hold on;
end
plot(x,YYY(1,:),'--','LineWidth',2);hold on;
xlabel(name,'Fontsize',20);
H=legend('GA | RR','GR | RA','GA | RA','GR | RR','\eta_0');
set(H,'Fontsize',18);
h=suptitle(['Fraction of cells under different ',name]);
set(h,'Fontsize',25);

if frac_type == 'GY'
    ylabel('Fraction of G+Y cells','Fontsize',20);
elseif frac_type == 'Y'
    ylabel('Fraction of Y cells','Fontsize',20);
end

end


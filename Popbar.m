function f=Popbar(x,name,frac_type,y0,y1,y2,y3,y4)

% name can be: [Kan],[Cm],Dilution rate/Hr
n=length(x);
for i=1:n
    Y(i,:)=[y0(i),y1(i),y2(i),y3(i),y4(i)];
end
figure;
b=bar(x,Y,1,'LineWidth',2);
ylim([0,1.0]);
set(gca,'LineWidth',2,'Fontsize',18);
xlabel(name,'Fontsize',20);
title('Ranking of 5 Scenarios');
legend({'\eta_0','GA|RR','GR|RA','GA|RA','GR|RR'},'Fontsize',18);
if frac_type == 'GY'
    ylabel('Fraction of Plasmid-carrying Cells','Fontsize',20);
elseif frac_type == 'Y'
    ylabel('Fraction of Transconjugants','Fontsize',20);
end
end


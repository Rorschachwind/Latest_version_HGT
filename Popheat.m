function f = Popheat(x,z,namex,namez,frac_type,y0,y1,y2,y3,y4)
%x is the variable in the second loop, z is the variable in the first loop
%plot imagesc for different combination
figure;
subplot(2,2,1);
imagesc(x,z,y1);
set(gca,'LineWidth',2,'Fontsize',18);
xlabel(namex,'Fontsize',20);
ylabel(namez,'Fontsize',20);
title('GA | RR','Fontsize',20);

subplot(2,2,2);
imagesc(x,z,y2);
set(gca,'LineWidth',2,'Fontsize',18);
title('GR | RA','Fontsize',20);

subplot(2,2,3);
imagesc(x,z,y3);
set(gca,'LineWidth',2,'Fontsize',18);
title('GA | RA','Fontsize',20);

subplot(2,2,4);
imagesc(x,z,y4);
set(gca,'LineWidth',2,'Fontsize',18);
title('GR | RR','Fontsize',20);
colorbar;

if frac_type == 'GY'
    h=suptitle('Franction of Plasmid-carrying Cells');
elseif grac_type == 'Y'
     h=suptitle('Franction of Plasmid-carrying Cells');
end
set(h,'Fontsize',25);
end
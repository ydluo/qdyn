clc;
ww=-p.Z(1:p.NX:end)./sin(p.DIP(1:p.NX:end)/180*pi)/1000;
xx=p.X(1:p.NX)/1000;
max_v_t=log10(max(ox1.v));
ifm=[50 52 60 66 68 74];
h1=figure(1)

subplot(4,2,[7:8])
for i=1:numel(ifm)
hold on;
  plot(ot1.t/year,log10(ot1.v),'b');
 % plot(ox1.t/year,mean_v_t,'g');
  plot(ox1.t(ifm(i))/year,max_v_t(ifm(i)),'ko');
 % plot(ox1.t(i)/year,mean_v_t(i),'go');
  legend('v_m_a_x');
  ylabel('log_1_0(V)');
  xlabel('Time (years)');
 % ylim([vvmin vvmax]);
  xlim([ox1.t(ifm(1))/year-0.3 ox1.t(ifm(end))/year+0.3]);
  hold off;
end

for i=1:6
 subplot(4,2,i)
 contourf(xx,ww,log10(reshape(ox1.v(:,ifm(i)),p.NX,p.NW))',40,'LineStyle','none');
 view(2)
 axis ij;
 caxis([-10 -1]);
 title(['t = ', num2str(ox1.t(ifm(i))/year,'%8.4f'),' year']);
 if i == 1 || i == 3 || i == 5
     ylabel('Along-dip (km)')
     
 end
 
 if i == 1 || i == 2 ||i == 3 || i ==4

     set(gca,'XTickLabel','')
 end
  if i == 6
     colorbar('South')
 end
 if i == 2 || i == 4
     set(gca,'XTickLabel','','YTickLabel','')
 end
 if i == 6
     xlabel('Along-strike (km)')
     set(gca,'YTickLabel','') 
 end
 
  if i == 5
     xlabel('Along-strike (km)')

 end
 
 axis equal
end

% for i=4:6
%  subplot(8,9,[28:30 37:39 46:48]+(i-4)*3)
%  contourf(xx,ww,log10(reshape(ox1.v(:,ifm(i)),p.NX,p.NW))',40,'LineStyle','none');
%  view(2)
%  axis ij;
%  caxis([-10 -1]);
%  title(['t = ', num2str(ox1.t(ifm(i))/year,'%8.4f'),' year']);
%  if i == 5
%     colorbar('South')
%     xlabel('Along-strike (km)')
%     set(gca,'YTickLabel','')
%  end
%   if i == 6
%     set(gca,'YTickLabel','')
%  end
% end



  print(h1,'-dpdf',['JP_3d_600_nn.pdf']);
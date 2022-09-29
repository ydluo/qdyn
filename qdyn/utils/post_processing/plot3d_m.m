clc;
name='JP_3d_L300_lc_rnd';
mname=[name '.avi'];
year=3600*24*365;
vvmin=-2;  %min v/vpl
vvmax=8;   %max v/vpl
t_const=1;  %time window for plotting v, if t_const=1, constant twin window, else, adjustable
t_win=10;  %time window for plotting v if t_const=1;
i_view=2;  %3 for 3d viw, 2 for 2d view
max_v_t=log10(max(ox1.v)./p.V_SS);
mean_v_t=log10(mean(ox1.v)./p.V_SS);
iplot=1;
for i=1:iplot:numel(ox1.t)
% for i=1:10
  disp(['Plotting ', num2str(i),'/',num2str(numel(ox1.t))]);
  h1=figure(1);
  subplot(6,6,[7:6:25]);
  temp_v=ox1.v(:,i);
  temp_v=reshape(temp_v,p.NX,p.NW);
  hold on;
  plot(log10(max(temp_v)./p.V_SS),p.Z(1:p.NX:end)/1000,'b');
  plot(log10(mean(temp_v)./p.V_SS),p.Z(1:p.NX:end)/1000,'g');
  legend('v_m_a_x','v_m_e_a_n');
  xlabel('log_1_0(V/V_p_l)');
  ylabel('Depth : km');
  hold off;
  subplot(6,6,[31:36]);
  hold on;
  plot(ox1.t/year,max_v_t,'b');
  plot(ox1.t/year,mean_v_t,'g');
  plot(ox1.t(i)/year,max_v_t(i),'bo');
  plot(ox1.t(i)/year,mean_v_t(i),'go');
  legend('v_m_a_x','v_m_e_a_n');
  ylabel('log_1_0(V/V_p_l)');
  xlabel('Time: years');
  ylim([vvmin vvmax]);
  if t_const == 1
      xlim([ox1.t(i)/year-t_win ox1.t(i)/year+t_win]);
  else
     if i<numel(ox1.t)
        xlim([ox1.t(i)-(ox1.t(i+1)-ox1.t(i))*5 ox1.t(i)+(ox1.t(i+1)-ox1.t(i))*5]/year);
     else
        xlim([ox1.t(i)-(ox1.t(i)-ox1.t(i-1))*5 ox1.t(i)+(ox1.t(i)-ox1.t(i-1))*5]/year); 
     end
  end
  hold off;
  subplot(6,6,[2:6 8:12 14:18 20:24 26:30]);
  scatter3(p.X/1000,p.Y/1000,p.Z/1000,[],log10(ox1.v(:,i)/p.V_SS),'s','filled');
  caxis([-2 8]);
  colorbar('NorthOutside');
  if i_view == 3
      axis equal;
      view(i_view);
  else
      view(i_view);  
  end
  title(['time = ', num2str(ox1.t(i)/year,'%15.8f'),' year']);
  zlabel('Depth : km');
  xlabel('Location along-strike: km');
  ylabel('Location Y: km');
%  print(h1,'-djpeg','-r1000',[name, num2str(i), '.jpg']);
  print(h1,'-dpdf',[name,'_m', num2str(i), '.pdf']);
%  mov(i)=getframe;
  clf(h1);
end
%movie2avi(mov,mname);

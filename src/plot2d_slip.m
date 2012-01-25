clc;
name='JP_2d_slip_z';
mname=[name '.avi'];
year=3600*24*365;
s0=215;
s1=874;
%for i=1:numel(ox1.t)
%for i=1:10
%  disp(['Plotting ', num2str(i),'/',num2str(numel(ox1.t))]);
  h1=figure(1);
  plot((ox1.d(:,s1)-ox1.d(:,s0)),p.Z/1000);
%  caxis([-2 8]);
%  colorbar('NorthOutside');
%  axis equal;
  title([ num2str(ox1.t(s0)/year,'%15.8f'),' to ', num2str(ox1.t(s1)/year,'%15.8f'),' year']);
  ylabel('Depth : km');
  xlabel('Slip: m');
%  ylabel('Location Y: km');
%  print(h1,'-djpeg','-r1000',[name, num2str(i), '.jpg']);
  print(h1,'-dpdf',[name, '.pdf']);
%  mov(i)=getframe;
  clf(h1);
%end
%movie2avi(mov,mname);


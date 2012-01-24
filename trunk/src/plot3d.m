
clc;
name='JP2D';
mname=[name '.avi'];
year=3600*24*365;
for i=1:numel(ox1.t)
%for i=1:10
  disp(['Plotting ', num2str(i),'/',num2str(numel(ox1.t))]);
  h1=figure(1);
  scatter3(p.X/1000,p.Y/1000,p.Z/1000,[],log10(ox1.v(:,i)/p.V_SS),'s','filled');
  caxis([-2 8]);
  colorbar;
  axis equal;
  title(['time = ', num2str(ox1.t(i)/year),' year']);
  zlabel('Depth : km');
  print(h1,'-djpeg','-r600',[name, num2str(i), '.jpg']);
%  mov(i)=getframe;
  clf(h1);
end
%movie2avi(mov,mname);


clc;
clear;
clf;
name='tremor3D_asp_nn_s1_s40';

year=3600*24*365;
day=3600*24;
hour=3600;
iplot=1;        %frame interval for output
%plot control
istart=1001;

icaxis_const=1; % 1 for constant caxis of slip rate   
icaxis_const_d=1; % 1 for constant caxis of slip  
iaxis_equal=1;	% 1 for equal axis  
% vvmin=0;  %min (v)
% vvmax=2;   %max (v)
% ddmin=0;  %min slip
% ddmax=60;   %max slip
i_view=2;  %3 for 3d viw, 2 for 2d view(fault view)
fps=12;  %framerate for movie, FPS
ppw=800;    %figure width
pph=300;    %figure height
t_win=10;


p = Qdyn_read_in();
plot_nw=floor(p.NW/2);

ss = Qdyn_read_seism('fort.22');

ox0=Qdyn_read_ox_seq(['fort.',num2str(istart)]);
p.X=ox0.X;
p.Y=ox0.Y;
p.Z=ox0.Z;
p.N=p.NW*p.NX;

h1=figure(1);
set(h1,'Position',[100 100 100+ppw 100+pph])

for iw=1:p.NW    
    cod(1+(iw-1)*p.NX:iw*p.NX)=p.Y(1+(iw-1)*p.NX:iw*p.NX)./cosd(p.DIP(iw));
end


h1=figure(1);
set(h1,'Position',[100 100 100+ppw 100+pph])
ylim([min(p.X)/1000 max(p.X)/1000]);
scatter(ss.t/day,p.X(ss.iloc)/1000,3*log10(ss.v/1e-5),cod(ss.iloc)/max(cod),'d','filled');
colormap('hot');
colorbar('west');
% set(p,'Markersize',log10(ss.v/1e-9))
xlabel('Time: (Days)','FontSize',16);
ylabel('Location along-strike: (km)','FontSize',16);
set(gca,'FontSize',16);
print(h1,'-depsc2',[name,'_proj_seis_strike', '.eps']);

h2=figure(2);
set(h2,'Position',[100 100 100+ppw 100+pph])
ylim([min(cod)/1000 max(cod)/1000]);
colormap('default');
scatter(ss.t/hour,cod(ss.iloc)/1000,3*log10(ss.v/1e-5),p.X(ss.iloc)/max(p.X),'d','filled');
colorbar('west');
xlabel('Time: (Hours)','FontSize',16);
ylabel('Location along-dip: (km)','FontSize',16);
set(gca,'FontSize',16);
print(h2,'-depsc2',[name,'_proj_seis_dip', '.eps']);


% %  mov(i)=getframe;
%   framev(i)=getframe(h1);
%   clf(h1);
% end
% writegif('V.gif',framev,1/fps);
% writegif('V_d.gif',framev(1:ceil(numel(framev)/120):end),1/24);

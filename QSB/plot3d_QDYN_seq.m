clc;
clear;
name='Dc005to007';
name2='test';

year = 365*24*3600;

iplot = [1001:10:1505];
% iplot = 20001;

iskip=1;       % # ponits skipped while plotting

icaxis_const=1; % 1 for constant caxis of slip rate  
vmin=1e-15;
vmax=0.1;        %m/s
icaxis_const_d=1; % 1 for constant caxis of slip
dmin=0;
dmax=8;        %m
icaxis_d_auto=0;    % =1 for auto setting dmax

i_gif = 1; % =1 to plot gif

iend_auto=0;        %=1 for auto end plotting while vmax< vth_plot;
vth_plot=0.001;     %m/s

iaxis_equal=1;	% 1 for equal axis  
% vvmin=0;  %min (v)
% vvmax=2;   %max (v)
% ddmin=0;  %min slip
% ddmax=60;   %max slip
fps=12;  %framerate for movie, FPS
ppw=1200;    %figure width
pph=200;    %figure height

isnap=max(iplot);
display(['Processsing Snapshot # ', num2str(isnap),'...']);
d = Qdyn_read_ox_seq(['fort.' num2str(isnap)]);

if icaxis_d_auto == 1
    dmax=max(d.d);
end 


h1=figure(1);
set(h1,'Position',[100 100 ppw pph])
set(h1,'Color',[1 1 1]);
scatter(d.X(1:iskip:end)/1000,d.Z(1:iskip:end)/1000,[],d.d(1:iskip:end),'s','filled');colorbar;
if iaxis_equal == 1
    axis equal;
end
axis tight;
if icaxis_const_d == 1
    caxis([dmin dmax]);
end
set(gca,'FontSize',16);
title('Slip (m)','FontSize',16);
xlabel('Along-strike: (km)','FontSize',16);
ylabel('Along-dip: (km)','FontSize',16);
text(max(d.X)/1000,min(d.Z)/1000,['time = ', num2str(d.t/year,'%15.5f'),' Years'],...
      'color','White','HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16); 
clf(h1); 
  
h2=figure(2);
set(h2,'Position',[100 300+pph ppw pph])
set(h2,'Color',[1 1 1]);
scatter(d.X(1:iskip:end)/1000,d.Z(1:iskip:end)/1000,[],log10(d.v(1:iskip:end)),'s','filled');colorbar;
if iaxis_equal == 1
    axis equal;
end
axis tight;
if icaxis_const == 1
    caxis([log10(vmin) log10(vmax)]);
end
set(gca,'FontSize',16);
title('Log10 slip rate (m/s)','FontSize',16);
xlabel('Along-strike: (km)','FontSize',16);
ylabel('Along-dip: (km)','FontSize',16);
text(max(d.X)/1000,min(d.Z)/1000,['time = ', num2str(d.t/year,'%15.5f'),' Years'],...
    'color','White','HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16); 
clf(h2);  
 

ifm=0;


for iii = 1:1:numel(iplot)
     
    isnap = iplot(iii);

    display(['Processsing Snapshot # ', num2str(isnap),'...']);
    d = Qdyn_read_ox_seq(['fort.' num2str(isnap)]);
   
    if (iend_auto ==1) && (max(d.v) <= vth_plot) && (isnap > 10)
        disp(['Event end at t = ', num2str(d.t/year,'%15.5f'),' Years', '  Snapshot # ', num2str(isnap)]);
        break
    end
    
    ifm=ifm+1;  

    h1=figure(1);
    set(h1,'Position',[100 100 ppw pph])
    set(h1,'Color',[1 1 1]);
    scatter(d.X(1:iskip:end)/1000,d.Z(1:iskip:end)/1000,[],d.d(1:iskip:end),'s','filled');colorbar;
    if iaxis_equal == 1
        axis equal;
    end
    axis tight;
    if icaxis_const_d == 1
        caxis([dmin dmax]);
    end
    set(gca,'FontSize',16);
    title('Slip (m)','FontSize',16);
    xlabel('Along-strike: (km)','FontSize',16);
    ylabel('Along-dip: (km)','FontSize',16);
    text(max(d.X)/1000,min(d.Z)/1000,['time = ', num2str(d.t/year,'%15.5f'),' Years'],...
      'color','White','HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16); 
    fD=getframe(h1);
    imwrite(fD.cdata,['D_',name,name2,'_',num2str(isnap), '.jpg'],'jpg');    
    frameD(ifm)=fD;
    clf(h1); 

    h2=figure(2);
    set(h2,'Position',[100 300+pph ppw pph])
    set(h2,'Color',[1 1 1]);
    scatter(d.X(1:iskip:end)/1000,d.Z(1:iskip:end)/1000,[],log10(d.v(1:iskip:end)),'s','filled');colorbar;
    if iaxis_equal == 1
        axis equal;
    end
    axis tight;
    if icaxis_const == 1
        caxis([log10(vmin) log10(vmax)]);
    end
    set(gca,'FontSize',16);
    title('Log10 slip rate (m/s)','FontSize',16);
    xlabel('Along-strike: (km)','FontSize',16);
    ylabel('Along-dip: (km)','FontSize',16);
    text(max(d.X)/1000,min(d.Z)/1000,['time = ', num2str(d.t/year,'%15.5f'),' Years'],...
        'color','White','HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16); 
%    print(h2,'-djpeg',['V_',name,name2,'_',num2str(isnap), '.jpg']);
%    print(h2,'-depsc2', ['V_',name,name2,'_',num2str(isnap), '.eps']);
    fV=getframe(h2);
    imwrite(fV.cdata,['V_',name,name2,'_',num2str(isnap), '.jpg'],'jpg');    
    frameV(ifm)=fV;
    clf(h2); 

end

if i_gif ==1
    writegif('V.gif',frameV,1/fps);
    writegif('V_d.gif',frameV(1:ceil(numel(frameV)/240):end),1/12);
    writegif('D.gif',frameD,1/fps);
    writegif('D_d.gif',frameD(1:ceil(numel(frameD)/240):end),1/12);
end

% Generate figures for January 2012 report
%
% User can set DataFile and MakeFigs

%DefaultDataFile = '/home/luoyd/qdyn_mm4/JP_2D_s_twm10000ts1000L1033.3914nx1W206.6783nw104dip14.mat';
DefaultDataFile = 'JP_2D_s_twm10000ts1000L1033.3914nx1W206.6783nw104dip14.mat';
DefaultMakeFigs = [1 1];

L=200e3;  % assumed along-strike rupture length (200 km)

if ~exist('DataFile'), DataFile=DefaultDataFile; end
if ~exist('MakeFigs'), MakeFigs = DefaultMakeFigs; end  

if ~exist('OldDataFile','var'), OldDataFile=DataFile; end
if ~strcmp(OldDataFile,DataFile) | ~exist('ox1','var') 
  disp(['Loading ' DataFile])
  load(DataFile); 
end
OldDataFile = DataFile;

%--- figure 1: 2D earthquake cycles
if MakeFigs(1)
  figure(1)

  subplot(4,4,[1 9])
  plot( (p.B-p.A)*100, p.Z/1e3,p.SIGMA/100e6, p.Z/1e3) % set double axis
  xlabel('(b-a) and \sigma')
  ylabel('Depth (km)') 
  legend('b-a (%)','\sigma / 10^8', 'Location','SE')
  grid on

  subplot(4,4,[2 11])
  plot( ox0.d(:,1:3:end), p.Z/1e3,'b')
  xlabel('Slip (m)')
  xlim([0 850])
  set(gca,'YTickLabel','') 
  title('Tohoku earthquake cycles (QDYN 2D)')

  subplot(4,4,[4 12])
  k=find(max(ox1.v)>1e-3);
  d0 = ox1.d(:,k(1));
  k=k(1:15:end);
  %plot( ox1.d(:,k(end))-ox1.d(:,k(1)), p.Z/1e3 )
  plot( ox1.d(:,k)-repmat(d0,1,length(k)), p.Z/1e3, 'b' )
  xlim([0 75])
  xlabel('Slip (m)')
  set(gca,'YAxisLocation','right')
  ylabel('Depth (km)') 
  potency = sum((ox1.d(:,k(end))-d0).*p.DW(:)) *L;
  Mw = (log10(p.MU*potency)-9.1)/1.5;
  title(['M_w ',sprintf('%0.1f',Mw),' megathrust earthquake'])
   
  subplot(4,4,[13 16])
  year = 365*24*3600;
  semilogy(ot0.t/year, ot0.pdot/dx*dw0*L *p.MU)
  xlim([0 10e3])
  ylim([1e10 1e22])
  xlabel('Time (years)')
  ylabel('Moment rate (N.m/s)')

end

if MakeFigs(2)
end

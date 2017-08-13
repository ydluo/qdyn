% a brittle asperity (velocity-weakening) surrounded by creep (velocity-strengthening)

Lasp = 7; % asperity "size" normalized by Lc = mu*Dc/((b-a)*sigma)
L = 5;    % total fault size normalized by Lasp
AB_RATIO = 0.9; % a/b at the center of the asperity. Outside: a/b = 2-AB_RATIO
RUN_OR_READ = 1; % 1 = run the simulation, 0 = read a pre-computed simulation 
RESOLUTION = 7;	% minimum number of nodes per Lb-length, dx/Lb>RESOLUTION, usually 9. 
  		% For large a/b and L/Lc, RESOLUTION=7 might be ok

%-----------

filename = 'example_1.mat';

year = 3600*24*365;
%if exist('qdyn')~=2, addpath ~/2D_RUPTURE/RATE_AND_STATE/qdyn/ ; end
p = qdyn('set');

% characteristic half-lengths
Lb = p.MU*p.DC/p.SIGMA/p.B;
Lnuc = 1.3774*Lb;
Lc = Lb/(1-AB_RATIO);
Linf = 1/pi *(1-AB_RATIO)^2 *Lb;

Lasp = Lasp*Lc;
L = L*Lasp;

if RUN_OR_READ

  p.L = L;
  %p.W = 1000*p.L;
  p.FINITE=0;
  p.N = 2^nextpow2(RESOLUTION*p.L/Lb); 
  dx=p.L/p.N;
  Lb_over_dx = Lb/dx
  p.ACC = 1e-10;
  
  p = qdyn('set',p);
  
  p.A = p.B *( 1 +(1-AB_RATIO)*(1-2*exp(-(p.X/Lasp*2).^6)) ); 
  
  p.TMAX = 8* year; % 60?
  p.NTOUT=100;
  p.NXOUT=1;
  p.NSTOP=0;
  %p.V_0 = 1.01*p.V_SS ;
  p.V_0 = ( 1.+0.01*exp( -(p.X/Lasp*2).^6 ) )*p.V_SS ;
  p.V_0 = p.V_0/mean(p.V_0)*p.V_SS;
  [p,ot1,ox1] = qdyn('run',p) ;
  
  
%  p.TMAX = 0.4*year; 
%  p.NTOUT=10;
  p.TMAX = 1.2*year;  
  p.NTOUT=50;

  p.V_0 = ox1.v(:,end)';
  p.TH_0= ox1.th(:,end)';
  %p.V_0 =  (ox1.v(:,end)+ox1.v(end:-1:1,end))/2;
  %p.TH_0=  (ox1.th(:,end)+ox1.th(end:-1:1,end))/2;
  [p,ot,ox]=qdyn('run',p);
  
  save(filename)

else
  load(filename)
end

%--- warmup cycles
if exist('ot1','var')
figure(1)
semilogy(ot1.t/year,ot1.v,ot1.t/year,ot1.pdot)
xlabel('Time (years)')
legend('V_{max}','PR_{max}','Location','Best')
title('Warming cycles')
end


%--- last cycle: snapshots
figure(2)
subplot(411)
semilogy(ot.t/year,ot.v/p.V_SS, ox.t/year,max(ox.v)/p.V_SS,'s')
xlabel('Time (years)')
ylabel('V_{max} / V_{pl}')

subplot(4,1,[2 3])
semilogy(ox.x/Lc,ox.v/p.V_SS)
ylabel('V / V_{pl}')
axis([[-1 1]*p.L/Lc/2 1e-2 1e3])

subplot(414)
if length(p.B)==1, p.B=repmat(p.B,1,p.N); end
plot(p.X/Lc,p.A, p.X/Lc,p.B,'--')
axis([[-1 1]*p.L/Lc/2 0.005 0.015])
legend('a','b')
xlabel('x / L_c')


%--- last cycle: contour plot
% NOTE: export with option -painters
figure(3)
subplot(4,4,[4:4:12])
semilogx(ot.v/p.V_SS,ot.t/year,ot.pdot/p.V_SS/p.L,ot.t/year)
axis([0.1 1e3 0 ot.t(end)/year])
set(gca,'XTick',[1 10 100],'XTickLabel',[1 10 100])
title('V / V_{pl}')
legend('V_{max}','<V>')

subplot(4,4,[ 2:4 6:8 10:12]-1)
contourf(ox.x/Lc,ox.t/year,log10(ox.v/p.V_SS)',40,'LineStyle','none');
ylabel('Time (years)')
title('log( V / V_{pl} )')
colorbar('North')

subplot(4,4,[13:15])
plot(p.X/Lc,p.A, p.X/Lc,p.B,'--')
set(gca,'XLim',[ -1 1]*p.L/Lc/2)
set(gca,'YLim',[ 0.005 0.015])
legend('a','b')
xlabel('x / L_c')


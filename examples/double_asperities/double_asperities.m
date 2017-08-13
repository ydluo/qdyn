% Interaction between two brittle asperities (velocity-weakening) surrounded by creep (velocity-strengthening)

Lasp = 4; % asperity "size" normalized by Lc = mu*Dc/((b-a)*sigma)
Lasp_2 = 2; % asperity "size" normalized by Lc = mu*Dc/((b-a)*sigma)
L = 5;    % total fault size normalized by Lasp
AB_RATIO = 0.6; % a/b at the center of the asperity. Outside: a/b = 2-AB_RATIO
RUN_OR_READ = 1; % 1 = run the simulation, 0 = read a pre-computed simulation 
RESOLUTION = 7;	% minimum number of nodes per Lb-length, dx/Lb>RESOLUTION, usually 9. 
  		% For large a/b and L/Lc, RESOLUTION=7 might be ok
THETA_LAW=1;

%-----------

filename = 'example_1.mat';

year = 3600*24*365;
%if exist('qdyn')~=2, addpath ~/2D_RUPTURE/RATE_AND_STATE/qdyn/ ; end
p = qdyn('set');
p.THETA_LAW=THETA_LAW;
% characteristic half-lengths
Lb = p.MU.*p.DC./p.SIGMA./p.B;
Lnuc = 1.3774*Lb;
Lc = Lb/(1-AB_RATIO);
Linf = 1/pi *(1-AB_RATIO)^2 *Lb;

Lasp = Lasp*Lc;
Lasp_2 = Lasp_2*Lc;

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
  
  b=p.B;
  a=p.A;
  
  p.B(1:1:p.N)=a/(2-AB_RATIO);
  p.A(1:1:p.N)=a;
  
  p.A = p.A - b*AB_RATIO*exp(-((p.X-(Lasp+Lasp_2)*0.6)/Lasp*2).^6); 
  p.A = p.A - b*AB_RATIO*exp(-((p.X+(Lasp+Lasp_2)*0.6)/Lasp_2*2).^6); 
  
  p.TMAX = 8 * year; % 60?
  p.NTOUT=1000;
  p.NXOUT=1;
  p.NSTOP=0;
  %p.V_0 = 1.01*p.V_SS ;
  p.V_0 = ( 1.+0.01*exp( -(p.X/Lasp*2).^6 ) )*p.V_SS ;
  p.V_0 = p.V_0/mean(p.V_0)*p.V_SS;
  [p,ot1,ox1] = qdyn('run',p) ;
  
  
%  p.TMAX = 0.4*year; 
%  p.NTOUT=10;
  p.TMAX = 2 *year;  
  p.NTOUT= 10;

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
semilogy(ot1.t/year,ot1.v)
xlabel('Time (years)')
legend('V_{max}','Location','Best')
title('Warm-up cycles')
end

figure(2)
subplot(4,4,[4:4:12])
semilogx(ot.v/p.V_SS,ot.t/year,ot.pdot/p.V_SS/p.L,ot.t/year)
axis([0.1 1e7 0 ot.t(end)/year])
set(gca,'XTick',[1  1e3 1e6],'XTickLabel',[1  1e3 1e6])
title('V / V_{pl}')
legend('V_{max}','<V>')

subplot(4,4,[ 2:4 6:8 10:12]-1)
contourf(ox.x,ox.t/year,log10(ox.v/p.V_SS)',40,'LineStyle','none');
ylabel('Time (years)')
%ylim([0.12 0.17])
title('log( V / V_{pl} )')
%caxis([-10 10])
colorbar('North')

subplot(4,4,[13:15])
plot(ox.x,p.A./p.B)
%set(gca,'XLim',[ -1 1]*p.L/Lc/2)
%set(gca,'YLim',[ 0.005 0.015])
axis tight;
ylim([min(p.A./p.B)*0.8 max(p.A./p.B)*1.2])
legend('a/b')
xlabel('x')

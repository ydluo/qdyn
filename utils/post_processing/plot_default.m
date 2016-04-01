h1=figure;
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
colorbar('North')

subplot(4,4,[13:15])
plot(ox.x,p.DC)
%set(gca,'XLim',[ -1 1]*p.L/Lc/2)
%set(gca,'YLim',[ 0.005 0.015])
axis tight;
legend('Dc')
xlabel('x')


print(h1,'-djpeg','-r600',[filename '.jpg']);

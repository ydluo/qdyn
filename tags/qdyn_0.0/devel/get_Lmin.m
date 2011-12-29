% Lmin = the minimum crack-length reached during nucleation 
function Lmin = get_Lmin(ot)
[vmax,imax]=max(ot.v); % before a dynamic event
knuc = find( ot.vdot(1:imax-1)>0 & diff(ot.cl(1:imax))<0 & ot.v(1:imax-1)<0.5*vmax);
kmin=knuc(end);
Lmin=ot.cl(kmin);

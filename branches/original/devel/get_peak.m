%---------
% get location of peak
% with quadratic interpolation
function [xp,sp] = get_peak(s,x)
[sm,im] = max(s);
if im>1 & im<length(s)
  s1 = 0.5*( s(im+1)-s(im-1) );
  s2 = s(im-1)-2*s(im)+s(im+1);
  dx = x(2)-x(1);
  xp = x(im)-s1/s2*dx;
  sp = s(im)-0.5*s1*s1/s2;
else
  xp = x(im);
  sp = s(im);
end

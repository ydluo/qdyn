% compute length from stress peaks
% assumes one stress peak on each side (+/-)
function cl = crack_length(stress,x)

[nx,nt]=size(stress);
cl = zeros(nt,1);
for it=1:nt,
  cl(it) = get_cl(stress(:,it),x);
end


%---------
function cl = get_cl(s,x)

% 1. assuming one peak on each side of a global minimum
%    --> does not work when two minima instead of a single centered minimum
%[sm,im]=min(s);
%xL = get_peak( s.*(x<=x(im)),x );
%xR = get_peak( s.*(x>=x(im)),x );

% 2. assuming one peak on each side of center (+/-)
n=length(x);
xL = get_peak( s(1:n/2),x(1:n/2) );
xR = get_peak( s(n/2:n),x(n/2:n) );

cl = xR-xL;


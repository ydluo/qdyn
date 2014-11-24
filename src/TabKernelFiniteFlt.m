% I = TabKernelFiniteFlt(nmax,name)
%
% PURPOSE:	Tabulate the correction term of the static kernel for a finite fault
%		For fault quasistatic and quasidynamic simulations 
%		this needs to be done only once, with a large nmax.
%		Then you can just read/load the output files.
%			
% INPUT:	nmax 	nb of tabulated values
%			should be >= nb of elements in L (real fault length)
%		name	output file name
%
% OUTPUT:	I(1:nmax) term in brackets in eq.40 of
%   			A.Cochard and J.R.Rice
%   			A spectral method for numerical elastodynamic fracture 
%   			analysis without spatial replication of the rupture event
%   			J.Mech.Phys.Solids, Vol.45, No.8, pp.1393-1418, 1997
%			I(n) = 2/pi* \int_0^{n \pi} sin(u)/u du
%			The static kernel is then:
%			K(k) = -mu/2 |k| I(|k|L/pi)
%			where L is the physical length of the fault
%			and k=pi*n/L is the wavenumber ON THE PADDED LENGTH 2L 
%
%			I(:) is also exported to an ascii file name.tab
%			and to a Matlab binary name.mat
%
% EXAMPLE  	TabKernelFiniteFlt(4096,'kernel_I_4096');
%
function I = TabKernelFiniteFlt(nmax,name)
atol = 1e-14;
I = zeros(nmax,1);
for n=1:nmax,
  I(n) = integral(@mysinc,(n-1)*pi,n*pi,'AbsTol',atol);  
    %discontinued using of deprecated function quadl, 
    %benchmarked with a differnce smaller than 9th significant digit 
    %Do we need a higher precision integral tool?
end
I = cumsum(I);
I = 2/pi*I(:);

save(strcat(name,'.tab'),'I','-ASCII','-DOUBLE');
save(strcat(name,'.mat'),'I');

%--
function s=mysinc(x)
iszero= (x==0);
s=sin(x)./(x+iszero) +iszero; % avoids division by zero

% qdyn_okada_kernel calculates Okada static stress kernels 
% The fault is decomposed into a mesh of rectangular elements
% The kernels contain the shear and normal stresses on all elements
% induced by unit slip on each element
%		
% [Ks,Kn] = qdyn_okada_kernel(NN,mu,lambda,X,Y,Z,DIP,DX,DW,fault_type)
% 
% INPUTS:
% NN		number of fault elements if the fault is generally non-planar
%		or [NX,NW]=[number of elements along strike, along dip]
%       	if the fault is planar or dip changes only with depth
% mu, lambda	elastic moduli
% X,Y,Z		coordinates of element centers
% 		For planar faults, elements must be ordered as (along-strike,along-dip)
%		otherwise no particular order is assumed
% DIP,DX,DW	dip angle, element size along strike, element size along dip. 
%		Size can be 1, NW or NN
% fault_type 	 1 = right-lateral strike-slip
% 		-1 = left-lateral strike-slip
% 		 2 = thrust dip-slip
%		-2 = normal dip-slip
%
% OUTPUTS:
% Ks,Kn		kernel matrices, size (NN,NN) or (NX*NW,NX*NW)
%
% NOTES:
% This function calls the Fortran code kernel.exe
% To compile kernel.exe, follow the instructions in kernel.f90
%
% EXAMPLES:
% See example_1.m

function [Ks,Kn] = qdyn_okada_kernel(NN,mu,lam,X,Y,Z,DIP,DX,DW,fault_type)

  mesh_mode = length(NN);
  if mesh_mode==1
    N = NN; 
  else
    NX = NN(1);
    NW = NN(2);
    N = NX*NW;
    if length(DIP)==1, DIP=ones(NW,1)*DIP; end
    if length(DX)==1, DX=ones(NX,1)*DX; end
    if length(DW)==1, DW=ones(NW,1)*DW; end
    DIP=repmat(DIP(:)',NX,1);
    DX=repmat(DX(:),1,NW);
    DW=repmat(DW(:)',NX,1);
  end

  if ~exist('fault_type','var'), fault_type=1; end
 
  fid=fopen('kernel.in','w');
  fprintf(fid,'%i\n' , fault_type);
  fprintf(fid,'%.15g %.15g\n', lam,mu);
  fprintf(fid,'%u\n' , mesh_mode);
  if mesh_mode==1
    fprintf(fid,'%u\n' , N);
  else
    fprintf(fid,'%u %u\n' , NX, NW);
  end
  fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g\n',...
        [X(:),Y(:),Z(:),DIP(:),DX(:),DW(:)]');
  fclose(fid);

  status = system('./kernel.exe');

  fid = fopen('kernel.s.out');
  Ks = fread(fid,[N,inf],'double');
  fclose(fid);
    
  fid = fopen('kernel.n.out');
  Kn = fread(fid,[N,inf],'double');
  fclose(fid);

  % Generate full kernel matrix
  if mesh_mode>1
    K0s = reshape(Ks, NX,NW,NW);
    K0n = reshape(Kn, NX,NW,NW);
    Ks = zeros(NX,NW,NX,NW);
    Kn = zeros(NX,NW,NX,NW);
    for jo=1:NW
    for io=1:NX
      for js=1:NW
      for is=1:NX
        Ks(is,js,io,jo) = K0s(1+abs(io-is),js,jo);
        Kn(is,js,io,jo) = K0n(1+abs(io-is),js,jo);
      end
      end
    end
    end
    Ks = reshape(Ks,N,N);
    Kn = reshape(Kn,N,N);
  end
 
end

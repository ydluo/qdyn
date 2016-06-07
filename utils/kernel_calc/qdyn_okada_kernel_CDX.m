% qdyn_okada_kernel calculates okada shear stress kernel 
%  for a mesh with constant strike and constant along-strike grid size dx
%		
function K = qdyn_okada_kernel_CDX(N,NW,NX,mu,lam,X,Y,Z,DIP,XX,WW)

  % Compute non-redundant elements of the kernel matrix
    fid=fopen('qdyn.in','w');
    fprintf(fid,'%u     meshdim\n' , 2);  
    fprintf(fid,'%u %u %u    N\n' , N,NW,NX);      
    fprintf(fid,'%.15g %.15g  LAM,MU\n', lam,mu);       
    fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g\n',...
        [X(:),Y(:),Z(:),DIP(:),XX(:),WW(:)]');
    fclose(fid);

    status = system('./kernel.exe');

    K0 = textread('kernel.out');
    
  % Generate full kernel matrix
    K = zeros(N);
    iiK = 0;
    % i:src,  j:obs
    for j= 1:N
        for i = 1:N
            iiK = iiK+1;
            isz = ceil(i/NX);
            isx = i - (isz-1)*NX;
            ioz = ceil(j/NX);
            iox = j - (ioz-1)*NX;
%            if mod(iiK,ceil(N*N/100)) == 0
%                disp([num2str(floor(iiK/N^2*100)) '%']);
%            end    
            II = N*(ioz-1) + NX*(isz-1) + 1 + abs(iox-isx);
            K(iiK) = K0(II);           
        end
    end

end

% qdyn_okada_kernel calculates okada shear stress kernel 
%  for a general fault mesh
%		
function [K] = qdyn_okada_kernel(N,mu,lam,X,Y,Z,DIP,XX,WW)

    fid=fopen('qdyn.in','w');
    fprintf(fid,'%u     meshdim\n' , 1);  
    fprintf(fid,'%u     N\n' , N);      
    fprintf(fid,'%.15g %.15g  LAM,MU\n', lam,mu);       
    fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g\n',...
        [X(:),Y(:),Z(:),DIP(:),XX(:),WW(:)]');
    fclose(fid);

    status = system('./kernel.exe');

    zz = textread('kernel.out');
    K = reshape(zz,N,N);
    
end

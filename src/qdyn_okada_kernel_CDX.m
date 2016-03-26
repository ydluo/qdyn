% qdyn_okada_kernel qdyn module to calculate okada stress kernel at giving
%     points
%		


function [K0] = qdyn_okada_kernel_CDX(N,NW,NX,mu,lam,X,Y,Z,DIP,XX,WW)


    fid=fopen('qdyn.in','w');
    fprintf(fid,'%u     meshdim\n' , 68);  
    fprintf(fid,'%u %u %u    N\n' , N,NW,NX);      
    fprintf(fid,'%.15g %.15g  LAM,MU\n', lam,mu);       
    fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g\n',...
        [X(:),Y(:),Z(:),DIP(:),XX(:),WW(:)]');
 
    fclose(fid);

    % solve
%     status = system('~/qdyn_svn/trunk/src/qdyn');
     status = system('./qdyn');
%    status = system('~/bin/qdyn');
    % rename input and output files

    K0 = textread('fort.68');
    
%    K = reshape(zzz,N,N);
    
end

% qdyn_okada_kernel qdyn module to calculate okada stress kernel at giving
%     points
%		


function [K] = qdyn_okada_kernel(N,mu,lam,X,Y,Z,DIP,XX,WW)


    fid=fopen('qdyn.in','w');
    fprintf(fid,'%u     meshdim\n' , 66);  
    fprintf(fid,'%u     N\n' , N);      
    fprintf(fid,'%.15g %.15g  LAM,MU\n', lam,mu);       
    fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g\n',...
        [X(:),Y(:),Z(:),DIP(:),XX(:),WW(:)]');
 
    fclose(fid);

    % solve
%     status = system('~/qdyn_svn/trunk/src/qdyn');
     status = system('./qdyn');
%    status = system('~/bin/qdyn');
    % rename input and output files

    zz = textread('fort.66');
    zzz = zz(:,1);
    K = reshape(zzz,N,N);
    
end

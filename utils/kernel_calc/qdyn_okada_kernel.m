% qdyn_okada_kernel calculates okada shear stress kernel 
%  for a general fault mesh
%		

function [K] = qdyn_okada_kernel(N,mu,lam,X,Y,Z,DIP,XX,WW)

    fid=fopen('kernel.in','w');
    fprintf(fid,'%u\n' , 1);  
    fprintf(fid,'%u\n' , N);      
    fprintf(fid,'%.15g %.15g\n', lam,mu);       
    fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g\n',...
        [X(:),Y(:),Z(:),DIP(:),XX(:),WW(:)]');
    fclose(fid);

    status = system('./kernel.exe');

    zz = textread('kernel.out');
    [xd,zd]=meshgrid(X,Z);
    K = reshape(zz(:,1),[N,N]);
    
end

%scatter3(X,Y,Z,[],K(200,:))

%zp=textread('kernel_py.in');
%Kp=reshape(zp(:,1),[N,N]);
%scatter3(X,Y,Z,[],Kp(200,:))
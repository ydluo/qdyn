clc;
clear;


x_off=-256e3;
y_off=0;
z_off=0;

isnap=60000;
%plot control
dt=0.005;
ft_id=1;
dir='./';
i_v_corr = 1;	% = 1 to set Vmin = v_corr
v_corr = 1e-20;

iskip=5;       % # ponits skipped while plotting


display(['Processsing Snapshot # ', num2str(isnap),'...']);
d = FSEM3D_snapshot(isnap,dir,ft_id);
d.X = d.X*1000;
d.Y = d.Y*1000;
d.Z = d.Z*1000;

if i_v_corr =  1
	d.Vx = max(d.Vx,v_corr);
end
disp('Done');

disp('Loading qdyn.in ...');
p=Qdyn_read_in();
disp('Done');

disp(['Loading QDYN output :  fort.20001 ...']);
o = Qdyn_read_ox_seq('fort.20001');
disp('Done');

p.X=o.X+x_off;
p.Y=o.Y+y_off;
p.Z=o.Z+z_off;

disp('Matching Values ...');

p.N = p.NX*p.NW;

for i=1:p.N
    
    [cc, id] = min((d.X(:)-p.X(i)).^2 + (d.Z(:)-p.Z(i)).^2);
    
    p.TH_0(i)=d.S(id);
    p.V_0(i)=d.Vx(id);


    if mod(i,ceil(p.N/100)) == 0
        disp([num2str(i/ceil(p.N/100)) '% Complete'])
    end
end

p.X=p.X-x_off;
p.Y=p.Y-y_off;
p.Z=p.Z-z_off;

disp('Done');

disp('Writing qdyn.in ...')
 % export qdyn.in
    fid=fopen('qdyn.in','w');
    fprintf(fid,'%u     meshdim\n' , p.MESHDIM); 
    if p.SIGMA_CPL == 1
      p.NEQS = 3;
    end
    if p.MESHDIM == 2 || p.MESHDIM == 4;
      fprintf(1, 'MESHDIM = %d\n', p.MESHDIM); %JPA should "1" be "fid" instead?
      fprintf(fid,'%u %u     NX, NW\n' , p.NX, p.NW);      
      fprintf(fid,'%.15g %.15g  %.15g      L, W, Z_CORNER\n', p.L, p.W, p.Z_CORNER);
      fprintf(fid,'%.15g %.15g \n', [p.DW(:), p.DIP_W(:)]');
    else  
      fprintf(fid,'%u     NN\n' , p.N);      
      fprintf(fid,'%.15g %.15g      L, W\n', p.L, p.W);
    end
    
    if p.MESHDIM == 1;
        fprintf(fid,'%u   finite\n', p.FINITE);
    end   
    
    fprintf(fid,'%u   itheta_law\n', p.THETA_LAW);
    fprintf(fid,'%u   i_rns_law\n', p.RNS_LAW);
    fprintf(fid,'%u   i_sigma_cpl\n', p.SIGMA_CPL);    
    fprintf(fid,'%u   n_equations\n', p.NEQS);
    fprintf(fid,'%u %u %u %u %u %u  ntout, nt_coord, nxout, nxout_DYN, ox_SEQ, ox_DYN\n', p.NTOUT,p.IC,p.NXOUT,p.NXOUT_DYN,p.OX_SEQ,p.OX_DYN);     
    fprintf(fid,'%.15g %.15g %.15g %.15g    beta, smu, lambda, v_th\n', p.VS, p.MU, p.LAM, p.V_TH);
    fprintf(fid,'%.15g %.15g    Tper, Aper\n',p.TPER,p.APER);
    fprintf(fid,'%.15g %.15g %.15g %.15g    dt_try, dtmax, tmax, accuracy\n',p.DTTRY,p.DTMAX,p.TMAX,p.ACC);
    fprintf(fid,'%u   nstop\n',p.NSTOP);
    fprintf(fid,'%u %u  DYN_FLAG, DYN_SKIP\n',p.DYN_FLAG,p.DYN_SKIP);
    fprintf(fid,'%.15g %.15g %.15g    M0, DYN_th_on, DYN_th_off\n', p.DYN_M,p.DYN_TH_ON,p.DYN_TH_OFF);

          
    fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n',...
        [p.SIGMA(:),p.V_0(:),p.TH_0(:),p.A(:),p.B(:),p.DC(:),p.V1(:),p.V2(:),p.MU_SS(:),p.V_SS(:),p.IOT(:),p.IASP(:)]');
 
    fclose(fid);


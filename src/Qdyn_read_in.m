%--------Qdyn_read_in
%--------read qdyn.in
%

function d = Qdyn_read_in()

fid=fopen('qdyn.in');
rline=fgetl(fid); rdat = sscanf(rline,'%f');
d.MESHDIM=rdat(1);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
d.NX=rdat(1);
d.NW=rdat(2);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
d.L=rdat(1);
d.W=rdat(2);
d.Z_CORNER=rdat(3);

d.DW=zeros(d.NW,1);
d.DIP=zeros(d.NW,1);
for i=1:d.NW
    rline=fgetl(fid); rdat = sscanf(rline,'%f');
    d.DW(i)=rdat(1);
    d.DIP_W(i)=rdat(2);
end
rline=fgetl(fid); rdat = sscanf(rline,'%f');
d.THETA_LAW=rdat(1);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
d.RNS_LAW=rdat(1);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
d.NEQUS=rdat(1);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
d.NTOUT=rdat(1);
d.IC=rdat(2);
d.NXOUT=rdat(3);
d.OX_SEQ=rdat(4);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
d.VS=rdat(1);
d.MU=rdat(2);
d.LAM=rdat(3);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
d.TPER=rdat(1);
d.APER=rdat(2);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
d.DTTRY=rdat(1);
d.DTMAX=rdat(2);
d.TMAX=rdat(3);
d.ACC=rdat(4);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
d.NSTOP=rdat(1);
fclose(fid);


dd=importdata('qdyn.in',' ',11+d.NW);
d.SIGMA=dd.data(:,1);
d.V_0=dd.data(:,2);
d.TH_0=dd.data(:,3);
d.A=dd.data(:,4);
d.B=dd.data(:,5);
d.DC=dd.data(:,6);
d.V1=dd.data(:,7);
d.V2=dd.data(:,8);
d.MU_SS=dd.data(:,9);
d.V_SS=dd.data(:,10);



return
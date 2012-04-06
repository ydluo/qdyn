%--------Qdyn_read_ox_seq
%--------read sequential ox output file
%

function d = Qdyn_read_ox_seq(filename)

fid=fopen(filename);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
d.it=rdat(1);
d.ivmax=rdat(2);
d.n=rdat(3);
d.t=rdat(4);


dd=importdata(filename,' ',2);
d.X=dd.data(:,1);
d.Y=dd.data(:,2);
d.Z=dd.data(:,3);
d.time=dd.data(:,4);
d.v=dd.data(:,5);
d.th=dd.data(:,6);
d.vd=dd.data(:,7);
d.dtau=dd.data(:,8);
d.dtaud=dd.data(:,9);
d.d=dd.data(:,10);



return
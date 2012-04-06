%--------Qdyn_read_ox_seq
%--------read sequential ox output file
%

function d = Qdyn_read_ot(filename)

dd=importdata(filename,' ',6);
d.t=dd.data(:,1);
d.locl=dd.data(:,2);
d.cl=dd.data(:,3);
d.p=dd.data(:,4);
d.pdot=dd.data(:,5);
d.vc=dd.data(:,6);
d.thc=dd.data(:,7);
d.omc=dd.data(:,8);
d.tauc=dd.data(:,9);
d.dc=dd.data(:,10);
d.xm=dd.data(:,11);
d.v=dd.data(:,12);
d.th=dd.data(:,13);
d.om=dd.data(:,14);
d.tau=dd.data(:,15);
d.d=dd.data(:,16);


return
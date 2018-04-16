function d = QDYN_snapshot(isnap)

qdfile=sprintf('fort.%u',isnap);
dd = importdata(qdfile,' ', 2);


d.X=dd.data(:,1);
d.Y=dd.data(:,2);
d.Z=dd.data(:,3);
d.T=dd.data(:,4);
d.V=dd.data(:,5);
d.TH=dd.data(:,6);
d.DV_V=dd.data(:,7);
d.TAU=dd.data(:,8);
d.DTAU_DT=dd.data(:,9);
d.D=dd.data(:,10);


clear dd

return



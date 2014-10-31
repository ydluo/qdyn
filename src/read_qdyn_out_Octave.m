% for Octave (by Martijn van den Ende, m.p.a.vandenende@uu.nl, Oct 30 2014)
function [ot,ox] = read_qdyn_out_Octave(namet,namex)

  fileID = fopen(namet);
  result = textscan(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'Delimiter', ',', 'HeaderLines', 6);
  fclose(fileID);
  
  ot.t = result(1){1,1};
  ot.locl = result(3){1,1};
  ot.cl = result(4){1,1};
  ot.p = result(5){1,1};
  ot.pdot = result(6){1,1};
  ot.vc = result(7){1,1};
  ot.thc = result(8){1,1};
  ot.omc = result(9){1,1};
  ot.tauc = result(10){1,1};
  ot.dc = result(11){1,1};
  ot.xm = result(12){1,1};
  ot.v = result(13){1,1};
  ot.th = result(14){1,1};
  ot.om = result(15){1,1};
  ot.tau = result(16){1,1};
  ot.d = result(17){1,1};
  ot.sigma = result(18){1,1};
  
  
  % snapshots
  fid=fopen(namex);
  cosa = textscan(fid,'%f %f %f %f %f %f %f %f %f', 'commentstyle','shell', 'delimiter', ' ');
  fclose(fid);
  ox.x = cosa(1){1,1};
  ox.t = cosa(2){1,1};
  ox.v = cosa(3){1,1};
  ox.th= cosa(4){1,1};
  ox.vd= cosa(5){1,1};
  ox.dtau = cosa(6){1,1};
  ox.dtaud = cosa(7){1,1};
  ox.d = cosa(8){1,1};
  ox.sigma = cosa(9){1,1};



clc;
clear;

isnap=0;
ft_id=1;
dir='./';

fname = 'nodesonfault';

display(['Processsing Snapshot # ', num2str(isnap),'...']);
d = FSEM3D_snapshot(isnap,dir,ft_id);
d.X = d.X*1000;
d.Y = d.Y*1000;
d.Z = d.Z*1000;
disp('Done');

disp(['Writing ' fname ' ...']);

disp('Done');



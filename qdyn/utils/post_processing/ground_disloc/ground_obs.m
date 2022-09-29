%testname = 'okada_test.mat';
%testname = 'Run3D_cmb_loop14t1nt200SSAF_wei_HeteDW_twm2000acc5ts100L24nx1W19.68dsse0.8nw267DD19.68sigma_co2.mat';
testname = 'Run3D_loop9t200nt500SSAF_wei_HeteDW_twm2000acc5ts100L24nx1W19.68dsse0.8nw267DD19.68sigma_co2.mat';
strike_angle = 90;		%strike angle with clockwise frtom north (degrees)

%xout_min = 0;
%xout_inc = 2;
%xout_max = 24;
%yout_min = -4;
%yout_inc = 2;
%yout_max = +4;

x = [-3:2:27];      %os location in km
y = [-15:2:15];

inc_plot = 1;	%plot_increment


tic;
disp(['Loading input file: ' testname]);
load(testname);
tlap = toc;
disp('Done');
disp(['Time elapsed: ' num2str(tlap) 's'] );


disp(['Preparing for Okada: ...']);

oxo = ox3d;

vvo = oxo.v(:,1:inc_plot:end);
tto = oxo.t(1:inc_plot:end);
nno = numel(tto);

% arrange input of okada
dxo = p.L/p.NX*ones(p.N,1);
dwo = repmat(p.DW(:)',p.NX,1); dwo = dwo(:);
zedgeo = -(p.Z' - dwo);
dipo = p.DIP';
sao = strike_angle*ones(p.N,1);
eoo = p.X';
noo = p.Y';


ddipo = zeros(p.N,1);
dteno = zeros(p.N,1);

%noo = dteno;

%       m(1) = fault length along the strike direction (m)
%       m(2) = fault width in dip direction (m)
%       m(3) = if sin(d) > 0, m(3) = depth to lower edge of fault (m)
%              if sin(d) < 0, m(3) = depth to upper edge of fault (m)
%              where d = dip angle
%       m(4) = dip angle, from the horizontal (degrees) 
%       m(5) = strike, clockwise from N (degrees)
%       m(6) = East offset of midpoint of lower edge from origin (m)
%       m(7) = North offset of midpoint of lower edge from origin (m)
%       m(8) = strike slip (negative for right lateral, 0 if N/A)
%       m(9) = dip slip:    if sin(2d) < 0, then m(9) > 0 --> normal slip
%                                       and m(9) < 0 --> reverse slip
%                           if sin(2d) > 0, then m(9) > 0 --> reverse
%                                       and m(9) < 0 --> normal
%               where d = dip angle
%               m(9) 0 if dip-slip N/A
%       m(10) = tensile motion (positive for opening, 0 if N/A)
%
%   Note: See Okada (1985) for details on geometry parameterization.
%


% Arragne output cood
% obs_grid
%x = xout_min: xout_inc : xout_max;
%y = yout_min: yout_inc : yout_max;
z = 0;
[X,Y]=meshgrid(x, y);
XY = [X(:) Y(:)];
NXY = size(XY, 1);
XYZ = [XY'; ones(1,NXY)*z];
XYZ = XYZ*1000.0;         % m

mu = p.MU;             % Pa
poisson = p.LAM/2/(p.LAM + p.MU);

vxo = zeros(NXY,nno);
vyo = vxo;
vzo = vyo;

tlap = toc;
disp('Done');
disp(['Time elapsed: ' num2str(tlap) 's'] );

for iio = 1:1:nno
        disp(['Processing ' num2str(iio) ' of ' num2str(nno) ' snapshot']);
	m = [dxo dwo zedgeo dipo sao eoo noo -vvo(:,iio) ddipo dteno];
%        m = [dxo dwo zedgeo dipo sao eoo noo -ones(size(dxo))  ddipo dteno];
	[U, D, S, flag]=disloc3d(m', XYZ, mu, poisson);
	vxo(:,iio) = U(1,:);
        vyo(:,iio) = U(2,:);
        vzo(:,iio) = U(3,:);
	tlap = toc;
	disp('Done');
	disp(['Time elapsed: ' num2str(tlap) 's'] );
end


save(['Ground_' testname],'-v7.3');

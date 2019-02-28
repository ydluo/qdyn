% QDYN		Quasi-dynamic earthquake cycles on a fault embedded in a
%		homogeneous, linear, isotropic elastic medium.
%               This is a Matlab wrapper for the Fortran code qdyn.f90
%		Features:
% 		+ rate-and-state friction coefficient
%   			Mu(V,Theta) = Muss + a*ln(V/Vss) + b*ln(Theta/Thetass)
%           	  or with transitions between velocity-weakening and strengthening
%		  via cut-off velocities:
%   			Mu = Muss - a*ln(1+V1/V) + b*ln(1+Theta*Dc/V2)
%	  	+ evolution equations for the state variable: ageing law
%   			dTheta/dT = 1-V*Theta/Dc
%		  or slip law
%   			dTheta/dT = -V*Theta/Dc * log(V*Theta/Dc)
%                 or ageing law in the self-accelerating approximation:
%   			dTheta/dT = -V*Theta/Dc
%		+ spatially non-uniform friction parameters a, b, Dc, v1, v2
%		  and initial conditions v(0), theta(0), sigma(0)
%		+ quasistatic stress balance with radiation damping, no inertia, no elastodynamics
%		+ geometries: spring-block, straight linear fault in a 1D elastic medium,
%		  or 2D fault surface with depth-dependent dip angle in a 3D elastic half-space
%		+ boundary/loading conditions
%			. (in 2D) slip can be spatially periodic along-strike
%			  and the fault is loaded by steady displacement at a distance W from the fault
%			  (crustal plane model, mimics a finite seismogenic depth W)
%		  	. or (in 2D and 3D) the fault area governed by rate-and-state friction has a
%			  finite size and is loaded by steady sliding on the rest of the fault
%
% SYNTAX	[pars,ot,ox] = qdyn(mode,[parsin],['Property',Value,...])
%
% INPUTS 	mode	'set'	gives the default parameter structure (pars),
%				overrides by fields present in structure parsin
%				or by Property/Value pairs
%			'write' writes qdyn input file only
%			'run'	sets parameters and runs a simulation
%			'read' 	reads parameters and outputs from a previous simulation
%		parsin	parameter structure to override the default parameters
%		'Prop' 	property to be set, a fieldname of the parameter structure
%		Value   of the above property, overrides default and parsin
%
%		The parameters that can be set through 'parsin' or 'Prop/Value' pairs are listed below.
%		Their default values can be obtained by running:
%			pars = qdyn('set')
%
%		Parameters defining the geometry of the problem and loading:
%		FAULT_TYPE loading geometry of the fault
%			 1 = strike-slip (right-lateral)
%			-1 = strike-slip (left-lateral)
%			 2 = thrust
%			-2 = normal
%		MESHDIM	dimension of the problem:
%			0 = spring-block system
%			1 = 1D fault in a 2D elastic medium
%			2 = 2D fault in a 3D medium
%		MU 	shear modulus (Pa)
%		LAM 	elastic modulus LAMBDA for 3D simulations (Pa)
%		VS 	shear wave velocity (m/s). If VS=0 radiation damping is turned off
%		D	damage level = 1 - (damaged shear modulus) / (intact shear modulus)
%		H	if D>0, half-thickness of the fault damage zone
%			if D=0, half-thickness of elastic slab bisected by fault
%		L 	fault length if MESHDIM=1
%		    	stiffness is MU/L if MESHDIM=0
%		FINITE	boundary conditions if MESHDIM=1
%			0 = fault is infinitely long but slip is spatially periodic with period L,
%			    loaded by steady displacement at distance W from the fault
%			1 = fault is infinitely long but only a segment of length L has
%			    rate-and-state friction, the rest has steady slip. If you get the
%			    error message "finite kernel is too small", create a larger kernel file
%			    using the function TabKernelFiniteFlt.m, update the file name in
%			    subroutine init_kernel_2D of src/fault_stress.f90, and recompile
%			2 = same as 0 but slip is symmetric relative to first element
%			3 = same as 1 but slip is symmetric relative to first element
%		W  	distance between displacement loading and fault if MESHDIM=1 and FINITE=0
%		DIP_W	dipping angle (degree). If depthdependent, values must be given
%			from deeper to shallower depth.
%		Z_CORNER fault bottom depth (m, negative down)
%		SIGMA_CPL  normal stress coupling (only for dipping faults)
%			0 = disable
%			1 = enable
%		APER 	amplitude of additional time-dependent oscillatory shear stress loading (Pa)
%		TPER 	period of oscillatory loading (s)
%
%		Rate-and-state friction parameters:
%		A  	amplitude of direct effect
%		B  	amplitude of evolution effect
%		DC 	characteristic slip distance (m)
%		MU_SS	reference steady-state friction coefficient
%		V_SS	reference steady-state slip velocity (m/s)
%		V_PL	slip rate loading (m/s)
%		RNS_LAW	type of rate-and-state friction law:
%			0 = original
%			1 = with cut-off velocities
%		V1 	cut-off velocity of direct effect (m/s)
%		V2 	cut-off velocity of evolution effect (m/s), controls the transition
%			from weakening to strengtherning when a<b. V2 should be <= V1
%		CO	Cohesion (Pa)
%
%		THETA_LAW type of evolution law for the state variable:
%			0 = ageing in the "no-healing" approximation
%			1 = ageing law
%			2 = slip law
%
%		Initial conditions:
%		SIGMA = effective normal stress (Pa)
%		V_0 = initial slip velocity (m/s)
%		TH_0 = initial state (m/s)
%
%		Discretization and accuracy parameters:
%   SOLVER   ODE solver mode
%     1 = Bulirsch-Stoer
%     2 = Runge-Kutta-Fehlberg
%		N	number of fault elements if MESHDIM=1. Must be a power of 2.
%		NX	number of fault elements along-strike in 3D. Must be a power of 2 if FFT used along-strike.
%		NW 	number of fault elements along-dip in 3D. Must be a power of 2 if FFT used along-dip.
%		DW 	along-dip length (m) of each element along-dip, from deeper to shallower
%		TMAX 	total simulation time (s)
%		NSTOP 	stopping criterion:
%			0 = stop at t=TMAX
%			1 = stop end of slip localization phase
%			2 = stop at first slip rate peak
%			3 = stop at v > TMAX
%		DTTRY 	first trial timestep (s)
%		DTMAX	maximum timestep (0=unrestricted)
%		ACC	solver accuracy
%
%		Output control parameters:
%		OX_SEQ 	type of snapshot outputs
%			0 = all snapshots in a single output file (fort.19)
%			1 = one output file per snapshot (fort.1001, ...)
%		NXOUT 	spatial interval (in number of elements) for snapshot outputs
%		NTOUT 	temporal interval (number of time steps) for snapshot outputs
%		NTOUT_OT 	temporal interval (number of time steps) for time series outputs
%		OX_DYN	output specific snapshots of dynamic events defined by thresholds
%			on peak slip velocity DYN_TH_ON and DYN_TH_OFF (see below)
%			0 = disable
%			1 = enable outputs for event #i:
%				event start: fort.19998+3i
%				event end: fort.19999+3i
%				rupture time: fort.20000+3i
%		NXOUT_DYN spatial interval (in number of elements) for dynamic snapshot outputs
% 		DYN_TH_ON peak slip rate threshold defining the beginning of a dynamic event
%		DYN_TH_OFF peak slip rate threshold defining the end of a dynamic event
%		IC 	index of selected element for time series output (ot)
%		IOT 	Indices of elements for additional time series outputs.
%			Set IOT(i) = 1 to enable time series outputs at the i-th element.
%			By default IOT(i)=-0 and this output is not done.
%               	Each element has a separate output file named fort.xxxxx
%			where xxxxx is an index (different than i) that starts at 10001
%			and is incremented by 1 for each selected element. For instance if
%			IOT=[0 0 1 1], the output of elements i=3 and i=4 are in files
%			fort.10001 and fort.10002, respectively.
%		IASP	Flags for elements. This is for convenient identification purposes only,
%			it does not affect computation and outputs.
%
%		Parameters for integration with SPECFEM3D dynamic code:
%		DYN_FLAG integration with dynamic code
%			0 = disable
%			1 = enable: stop QDYN at the DYN_SKIP+1-th dynamic event
%			    with seismic moment > DYN_M
%		DYN_M	target seismic moment of a dynamic event
%		DYN_SKIP number of dynamic events to skip (warm up cycles)
%
%		Other parameters:
%		EXEC_PATH path to the Fortran qdyn executable.
%			The default path is the directory containing qdyn.m
%		NPROCS number of processors for parallel MPI runs
%			The default value is 1 (serial run)
%
% OUTPUTS 	pars	structure containing the parameters listed above, and:
%			X,Y,Z = coordinates of center of fault elements
%			DIP = dip angle of fault elements
%		ot	structure containing time series outputs
%			ot.t	output times
%			ot.locl	localization length (distance between stressing rate maxima)
%			ot.cl	crack length (distance between slip rate maxima)
%			ot.p	seismic potency
%			ot.pdot	seismic potency rate
%			-- outputs at the fault location with maximum slip rate:
%			ot.xm	location of maximum slip rate point
%			ot.v	maximum slip rate
%			ot.th	state variable theta
%			ot.om 	(slip rate)*theta/dc
%			ot.tau	stress
%			ot.d	slip
%			-- outputs at selected fault node with index pars.ic
%			ot.vc	slip rate
%			ot.thc	state variable theta
%			ot.omc 	slip_rate*theta/dc
%			ot.tauc	shear stress
%			ot.dc	slip
%		ox	structure containing snapshot outputs
%			ox.x	fault coordinates (for 2D)
%			ox.t 	output times
%			ox.v	slip rate
%			ox.th	state variable
%			ox.vd	slip acceleration
%			ox.dtau shear stress relative to initial value
%			ox.dtaud shear stress rate
%			ox.d 	slip
%			ox.sigma effective normal stress
%
% EXAMPLE	A run with initial velocity slightly above the default steady state:
%			p = qdyn('set');
%			[p,ot,ox] = qdyn('run','V_0',1.01*p.V_SS);
%			semilogy(ot.t,ot.v)
%
% AUTHORS	QDYN development team https://github.com/ydluo/qdyn#developers

function [pars,ot,ox] = qdyn(mode,varargin)

% NOTE ON VARIABLE NAMING CONVENTIONS
%	lower_case 	= local variables
%	UPPER_CASE 	= variables that will be wrapped into output structure 'pars'

%-- useful units
day = 60*60*24;
month = 30*day;
year = 365*day;

%--------- DEFAULT PARAMETERS ------------------------------------

NPROCS = 1; % Default serial, no-MPI.

MESHDIM=1;
NAME ='';	% title for the simulation

scriptName = mfilename('fullpath');
EXEC_PATH = fileparts(scriptName);	% default is same directory as qdyn.m

%-- fault and solver type
FAULT_TYPE = 1;
SOLVER = 1;

%-- medium
L= 2e3; 	% fault length (L scales the stiffness for the spring-block case)
FINITE=0;	% along strike bcs: 1=fixed-length asperity surrounded by steady creep, 0=periodic
W= 50e3;   	% out-of-plane dimension (ignored in spring-block)
MU = 30e9;	% shear modulus
LAM = 30e9;
VS = 3000; 	% shear wave velocity (if VS=0: turn off radiation damping)
D = 0;		% damage level
H = 0;		% half-thickness of the fault damage zone
V_TH= 1e-5; % threshold velocity for seismic events;

%-- numerical settings
N=1024/2; 	% number of grid cells
NX=100/2;
NW=10;
DW=1e3;
DIP_W=90.0;

Z_CORNER=-50e3;
IC=1;         % index of nodes for output ot
TMAX = 6*month;  % total simulation time
NSTOP = 0;	% stop at (0) tmax, (1) end of localization or (2) max slip rate
DTTRY = 1e2;   % first trial timestep
DTMAX = 0;	% maximum timestep (0=unrestricted)
ACC = 1e-7;     % solver accuracy
NXOUT = 8;	% space stride (cells) for snapshot outputs
NXOUT_DYN = 1;	% space stride (cells) for dynamic snapshot outputs
NTOUT = 100; 	% time stride (iterations) for snapshot outputs
NTOUT_OT = 1;   % time stride (iterations) for time series outputs
OX_SEQ = 0; 	% = 1 ; enable sequential ox output , from fort.1000 ...
OX_DYN = 0; % = 1 ; enable sequential snapshot of dynamic events
IOT = 0;    % = 1 to output ot of indicated nodes
IASP = 0;   % =1 indicating that it is an asperity
%-- friction
A = 0.9e-2;
B = 1e-2;
DC = 4e-4;
MU_SS = 0.6;
V_SS = 1e-9;
V_PL = V_SS;
THETA_LAW = 1;
RNS_LAW = 0;
SIGMA_CPL = 0;
CO = 0;
%-- initial conditions
SIGMA=1e8;
V_0=V_SS ;
TH_0=DC/V_SS;
V1=0.01;
V2=1e-7;
%-- branching fault
%JPA This is an undocumented feature implemented by Percy
BRANCH='.false.';
STRIKE=0;
XC=0; % Junction location.
YC=0; %

%-- periodic loading
APER = 0;
TPER = 1*year;

% dynamic intergrating parameters
DYN_FLAG = 0;
DYN_M = 1.e18; % M0= 10^(1.5Mw+9) [M0=1.e18 => Mw=6.0; M0=1.e21 => Mw=8.0]
DYN_SKIP = 0;
DYN_TH_ON = 1e-3;
DYN_TH_OFF = 1e-4;

%--------- INITIALIZE -------------------------------------

% override UPPER CASE VARIABLES with inputs
Parse_Inputs(varargin{:});

% generate the mesh
[X,Y,Z,DIP] = generate_mesh();

if MESHDIM<2 && NPROCS>1 
  disp('MPI parallelization is only implemented for MESHDIM=2. Resetting NPROCS=1.')
  NPROCS = 1;
end

if MESHDIM==1 & N < 2^nextpow2(N)
  error('N must be a power of 2')
end

% set steady state theta
TH_SS = DC./V_SS;

% wrap UPPER CASE VARIABLES in parameter structure fields with the same name
fpars = who;
for k= find( strcmp(fpars,upper(fpars)) )'
  pars.(fpars{k}) = eval(fpars{k}) ;
end

switch mode

case 'set'
  ot=[];
  ox=[];
  return % pars = qdyn('set',...)  --> do not compute, exit here

case 'read'
  pars = read_qdyn_in(NAME);
  pars.NAME = NAME;
  [ot,ox] = read_qdyn_out(NAME, pars.OX_SEQ);

case {'run', 'write'}

  % make vectors if constants
   DW(1:NW) = DW;
   DIP_W(1:NW) = DIP_W;
   A(1:N) = A;
   B(1:N) = B;
   DC(1:N) = DC;
   V_0(1:N) = V_0;
   TH_0(1:N)= TH_0;
   V1(1:N) = V1;
   V2(1:N) = V2;
   SIGMA(1:N) = SIGMA;
   MU_SS(1:N) = MU_SS;
   V_SS(1:N) = V_SS;
   V_PL(1:N) = V_PL;
   IOT(1:N) = IOT;
   IASP(1:N) = IASP;
   CO(1:N) = CO;

  % For branching faults.
  % NOTE: this is an undocumented feature implemented by Percy
   if strcmp(BRANCH,'.true.')
     export_branch_input();
     ot = [];
     ox = [];
     return;
   end

 % export list of recording points ("stations") on the fault
  if MESHDIM==2, export_stations(); end

 % export input file qdyn.in
  export_main_input()

  if strcmp(mode, 'write')
    ot = [];
    ox = [];
    return;
  end

 % solve
  if NPROCS==1
    status = system([EXEC_PATH filesep 'qdyn']);
  else
    status = system(['mpirun -np ' num2str(NPROCS) ' ' EXEC_PATH filesep 'qdyn']);
  end

 % rename input and output files
  if ~isempty(NAME)
    movefile('fort.18',[NAME '.ot']);
    movefile('fort.19',[NAME '.ox']);
    copyfile('qdyn.in',[NAME '.in']);
  end

 % output
  if NPROCS>1 % MPI parallel
    [ot,ox]= read_qdyn_out_mpi(NAME,OX_SEQ);
  else
    [ot,ox]= read_qdyn_out(NAME,OX_SEQ);
  end

otherwise
  error('mode must be: set, read or run')

end

%--
% This is a nested function: it has access to variables in the parent function
function [X,Y,Z,DIP]=generate_mesh()

  X = [];
  Y = [];
  Z = [];
  DIP = [];

  switch MESHDIM
  case 0
    X = [];
  case 1
    X = (-N/2+0.5:N/2-0.5) *L/N;
  case 2
   % set x, y, z, dip of first row
    cd = cos(DIP_W(1)/180.*pi);
    sd = sin(DIP_W(1)/180.*pi);
    X = 0. + (0.5+(0:NX-1))*L/NX;
    Y(1:NX) = 0.+0.5*DW(1)*cd;
    Z(1:NX) = Z_CORNER+0.5*DW(1)*sd;
    DIP(1:NX) = DIP_W(1);
   % set x, y, z, dip of row 2 to nw
    for i = 2:NW
      cd0 = cd;
      sd0 = sd;
      cd = cos(DIP_W(i)/180.*pi);
      sd = sin(DIP_W(i)/180.*pi);
      j0 = (i-1)*NX;
      X(j0+1:j0+NX) = X(1:NX);
      Y(j0+1:j0+NX) = Y(j0) + 0.5*DW(i-1)*cd0 + 0.5*DW(i)*cd;
      Z(j0+1:j0+NX) = Z(j0) + 0.5*DW(i-1)*sd0 + 0.5*DW(i)*sd;
      DIP(j0+1:j0+NX) = DIP_W(i);
    end
  otherwise
    error('MESHDIM must be 0, 1 or 2');
  end

end

%--
% This is a nested function: it has access to variables in the parent function
function export_branch_input()

  r=sqrt(X.^2+Y.^2);
  X=r.*cosd(STRIKE)+XC;
  Y=r.*sind(STRIKE)+YC;
  fid=fopen('qdyn_branch.in','w');
  fprintf(fid,'%d\n',N);
  fprintf(fid,'%15.6f\n',DIP_W(1));
  fprintf(fid,'%20.6f %20.6f\n',LAM,MU);
  fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n',...
    [X;Y;Z;SIGMA;V_0;TH_0;A;B;DC;V1;V2;MU_SS;V_SS;CO;V_PL]);
  fclose(fid);

end

%--
% This is a nested function: it has access to variables in the parent function
function export_stations()

   fid=fopen('stations.dat','w');
   nsta=1; %NOTE: For now one station but later will be extended to more stations
   fprintf(fid,'%d\n',nsta);
   fprintf(fid,'%.15g %.15g %.15g\n',X(IC),Y(IC),Z(IC));
   fclose(fid);

end

%--
% This is a nested function: it has access to variables in the parent function
function export_main_input()

  for iproc=1:NPROCS

    if NPROCS>1
      filename = ['qdyn' sprintf('%06i',iproc-1) '.in'];
     % MPI domain partitioning is done here:
      nw = floor(NW/NPROCS);
      i0 = (iproc-1)*NX*nw;
      iw0 = (iproc-1)*nw;
      % If NW/NPROCS is not integer, leave the rest to the last processor
      if iproc==NPROCS, nw = nw + mod(NW,NPROCS); end
      iloc = i0 + [1:NX*nw];
      iwloc = iw0 + [1:nw];
    else
      filename = 'qdyn.in';
      nw = NW;
      iloc = [1:N];
      iwloc = [1:NW];
    end

    fid = fopen(filename,'w');
    fprintf(fid,'%u     meshdim\n' , MESHDIM);
    if MESHDIM == 2
      fprintf(fid,'%u %u     NX, NW\n' , NX, nw);
      fprintf(fid,'%.15g %.15g %.15g      L, W, Z_CORNER\n', L, W, Z_CORNER);
      fprintf(fid,'%.15g %.15g \n', [DW(iwloc);DIP_W(iwloc)]);
    else
      fprintf(fid,'%u     NN\n' , N);
      fprintf(fid,'%.15g %.15g      L, W\n', L, W);
    end
    if MESHDIM == 1, fprintf(fid,'%u   finite\n', FINITE); end
    fprintf(fid,'%u   itheta_law\n', THETA_LAW);
    fprintf(fid,'%u   i_rns_law\n', RNS_LAW);
    fprintf(fid,'%u   i_sigma_cpl\n', SIGMA_CPL);
    fprintf(fid,'0 0 0     stress_coupling, thermal press., localisation\n');
    fprintf(fid,'%u %u %u %u %u %u %u  ntout_ot, ntout, nt_coord, nxout, nxout_DYN, ox_SEQ, ox_DYN\n', ...
                NTOUT_OT, NTOUT,IC,NXOUT,NXOUT_DYN,OX_SEQ,OX_DYN);
    fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g   beta, smu, lambda, D, H, v_th\n', ...
                VS, MU, LAM, D, H, V_TH);
    fprintf(fid,'%.15g %.15g    Tper, Aper\n',TPER,APER);
    fprintf(fid,'%.15g %.15g %.15g %.15g    dt_try, dtmax, tmax, accuracy\n',DTTRY,DTMAX,TMAX,ACC);
    fprintf(fid,'%u   nstop\n',NSTOP);
    fprintf(fid,'%u %u  DYN_FLAG, DYN_SKIP\n',DYN_FLAG,DYN_SKIP);
    fprintf(fid,'%.15g %.15g %.15g    M0, DYN_th_on, DYN_th_off\n', DYN_M,DYN_TH_ON,DYN_TH_OFF);
    fprintf(fid,'%u %u    FAULT_TYPE, SOLVER\n', FAULT_TYPE, SOLVER);

    fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %u %u %.15g %.15g\n',...
      [SIGMA(iloc);V_0(iloc);TH_0(iloc);A(iloc);B(iloc);DC(iloc);V1(iloc);V2(iloc); ...
       MU_SS(iloc);V_SS(iloc);IOT(iloc);IASP(iloc);CO(iloc);V_PL(iloc)]);

    if NPROCS>1
      fprintf(fid,'%.15g %.15g %.15g %.15g\n', ...
              [X(iloc);Y(iloc);Z(iloc);DIP(iloc)]);
    end

    fclose(fid);
  end
end % of nested function

end % of function qdyn

%-----------------------------------------------------------------------------
%
% PARSE_INPUTS sets variables in the caller function according to inputs.
%	By convention the variables to set have UPPER CASE names.
%
% WARNING: no checks yet
% Field names in parsin and Property's should match field names in pars
% or else the default value will be taken.
function Parse_Inputs(varargin)

if isempty(varargin)
    return
end

% 1. Override defaults by input structure
if isstruct(varargin{1})
 % process input parameter structure
  parsin = varargin{1};
  fparsin = fieldnames(parsin);
  for k=1:length(fparsin)
    assignin('caller', upper(fparsin{k}), parsin.(fparsin{k}) );
  end
 % separate the input Property/Value pairs
  varargin(1)=[];
end

% 2. Override defaults and previously set parameters by Property/Value pairs
for k = 2:2:length(varargin)
  assignin('caller', upper(varargin{k-1}), varargin{k} );
end
end % of function Parse_Inputs

%-----------
function pars = read_qdyn_in(name)
if ~exist('name','var') || isempty(name)
    name = 'qdyn';
end
name = [name '.in'];

if ~exist(name,'file') %
    name = [name(1:end-3) '000000.in'];  %MPI first node
    if exist(name,'file')
        disp(['ERROR: MPI input find: ' name '... cannot process']);
        return
    else
        disp('ERROR: no QDYN input');
        return
    end
end

fid = fopen(name);
MESHDIM = sscanf(fgetl(fid), '%u');
if MESHDIM == 2
    rdat = sscanf(fgetl(fid), '%u %u');
    NX = rdat(1); NW = rdat(2);
    rdat = sscanf(fgetl(fid), '%f %f %f');
    L = rdat(1); W = rdat(2); Z_CORNER = rdat(3);
    if ~exist('iwloc','var') || isempty(name)
	    iwloc = [1:NW];
    end
    rdat = textscan(fid, '%f %f'); DW(iwloc) = rdat(1); DIP_W(iwloc) = rdat(2);
else
    N = sscanf(fgetl(fid), '%u');
    rdat = sscanf(fgetl(fid), '%f %f'); L = rdat(1); W = rdat(2);
end
if MESHDIM == 1
    FINITE = sscanf(fgetl(fid), '%u');
end
THETA_LAW = sscanf(fgetl(fid), '%u');
RNS_LAW = sscanf(fgetl(fid), '%u');
SIGMA_CPL = sscanf(fgetl(fid), '%u');
rdat = sscanf(fgetl(fid), '%u %u %u %u');
rdat = sscanf(fgetl(fid), '%u %u %u %u %u %u %u'); NTOUT_OT = rdat(1); NTOUT = rdat(2); IC = rdat(3); NXOUT = rdat(4); NXOUT_DYN = rdat(5); OX_SEQ  = rdat(6); OX_DYN = rdat(7);
rdat = sscanf(fgetl(fid), '%f %f %f %f %f %f'); VS = rdat(1); MU = rdat(2); LAM = rdat(3); D = rdat(4); H = rdat(5); V_TH = rdat(6);
rdat = sscanf(fgetl(fid), '%f %f'); TPER = rdat(1); APER = rdat(2);
rdat = sscanf(fgetl(fid), '%f %f %f %f'); DTTRY = rdat(1); DTMAX = rdat(2); TMAX = rdat(3); ACC = rdat(4);
NSTOP = sscanf(fgetl(fid), '%u');
rdat = sscanf(fgetl(fid), '%u %u'); DYN_FLAG = rdat(1); DYN_SKIP = rdat(2);
rdat = sscanf(fgetl(fid), '%f %f %f'); DYN_M = rdat(1); DYN_TH_ON = rdat(2); DYN_TH_OFF = rdat(3);
rdat = sscanf(fgetl(fid), '%u %u'); FAULT_TYPE = rdat(1); SOLVER = rdat(2);
rdat = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %u %u %f %f');
SIGMA = rdat{1};
V_0 = rdat{2};
TH_0 = rdat{3};
A = rdat{4};
B = rdat{5};
DC = rdat{6};
V1 = rdat{7};
V2 = rdat{8};
MU_SS = rdat{9};
V_SS = rdat{10};
IOT = rdat{11};
IASP = rdat{12};
CO = rdat{13};
V_PL = rdat{14};
fclose(fid);

% wrap UPPER CASE variables in parameter structure fields with the same name
fpars = who;
for k = find( strcmp(fpars,upper(fpars)) )'
  pars.(fpars{k}) = eval(fpars{k});
end
end % of function read_qdyn_in

% read outputs from qdyn.f
function [ot,ox] = read_qdyn_out_mpi(name,OX_SEQv)

if exist('name','var') && ~isempty(name)
  namet = [name '.ot'];
  namex = [name '.ox'];
else
  namet = 'fort.18';
  namex = 'fort.19';
end

  % time series
  [ot.t,ot.vc, ot.thc, ot.omc, ot.tauc, ot.dc ] = ...
    textread(namet,'','headerlines',4);

  % snapshots
 if ~OX_SEQv 
  fid=fopen(namex);
  NSX=fscanf(fid,'# nx=%u');
  fclose(fid);
  cosa = textread(namex,'','commentstyle','shell');
  ncosa = size(cosa);
  NST=ncosa(1)/NSX;
  cosa=reshape(cosa,NSX,NST,ncosa(2));
  ox.x = cosa(:,1,1);
  ox.y = cosa(:,1,2);
  ox.z = cosa(:,1,3);

  ox.t = cosa(1,:,4)';
  ox.v = cosa(:,:,5);
  ox.th= cosa(:,:,6);
  ox.dtau = cosa(:,:,7);
  ox.dtaud = cosa(:,:,8);
  ox.d = cosa(:,:,9);
  ox.sigma = cosa(:,:,10);
 else 
     ox=[];
 end
end % of function read_qdyn_out_mpi

%-----------
% read outputs from qdyn.f
function [ot,ox] = read_qdyn_out(name,OX_SEQv)
    if exist('name','var') && ~isempty(name)
        namet = [name '.ot'];
        namex = [name '.ox'];
    else
        namet = 'fort.18';
        namex = 'fort.19';
    end

    if uimatlab
        % time series
        [~, file_sys] = system(['file ',namet]);
        isasc = textscan(file_sys, '%s', 2);
        if strcmp('ASCII', isasc{1}(2))
            [ot.t, ot.locl, ot.cl, ot.p, ot.pdot, ...
                ot.vc, ot.thc, ot.omc, ot.tauc, ot.dc, ...
                ot.xm, ot.v, ot.th, ot.om, ot.tau, ot.d, ot.sigma ] = ...
                textread(namet, '', 'headerlines', 6);
        else % binary output
            fid = fopen(namet, 'r');
            hdr = fread(fid, 'double');
            fclose(fid);

            hdr = reshape(hdr, 17, length(hdr)/17)';
            ot.t = hdr(:,1);
            ot.locl = hdr(:,2);
            ot.cl = hdr(:,3);
            ot.p = hdr(:,4);
            ot.pdot = hdr(:,5);
            ot.vc = hdr(:,6);
            ot.thc = hdr(:,7);
            ot.omc = hdr(:,8);
            ot.tauc = hdr(:,9);
            ot.dc = hdr(:,10);
            ot.xm = hdr(:,11);
            ot.v = hdr(:,12);
            ot.th = hdr(:,13);
            ot.om = hdr(:,14);
            ot.tau = hdr(:,15);
            ot.d = hdr(:,16);
            ot.sigma = hdr(:,17);
        end

        %snapshots
      %PG, Temporal. 
      ox=[];
      if ~OX_SEQv 
        [~, file_sys] = system(['file ',namex]);
        isasc = textscan(file_sys, '%s', 2);
        if strcmp('ASCII', isasc{1}(2))
            fid = fopen(namex);
            NSX = fscanf(fid, '# nx=%u');
            fclose(fid);
            cosa = textread(namex,'','commentstyle','shell');
            ncosa = size(cosa);
            NST = ncosa(1)/NSX;
            cosa = reshape(cosa, NSX, NST, ncosa(2));
            ox.x = cosa(:,1,1);
            ox.t = cosa(1,:,2)';
            ox.v = cosa(:,:,3);
            ox.th = cosa(:,:,4);
            ox.dtau = cosa(:,:,5);
            ox.dtaud = cosa(:,:,6);
            ox.d = cosa(:,:,7);
            ox.sigma = cosa(:,:,8);
        else % binary output
            fid = fopen(namex, 'r');
            nx = fread(fid, 1, 'int32');
            ox.x = fread(fid, nx, 'double');
            i = 1;
            while 1
                time = fread(fid, 1, 'double');
                if isempty(time)
                    fclose(fid);
                    break
                else
                    ox.t(i,:) = time;
                    hdr = fread(fid, [6, nx], 'double')';
                    ox.v(:,i) = hdr(:,1);
                    ox.th(:,i) = hdr(:,2);
                    ox.dtau(:,i) = hdr(:,3);
                    ox.dtaud(:,i) = hdr(:,4);
                    ox.d(:,i) = hdr(:,5);
                    ox.sigma(:,i) = hdr(:,6);
                    i = i + 1;
                end
            end
        end
      end	
    else
        [ot,ox] = read_qdyn_out_Octave(namet,namex)
    end
end % of function read_qdyn_out

%---
% adapted from http://www.mathworks.com/matlabcentral/fileexchange/23868-is-this-matlab-or-octave-
function uiIsMatLab = uimatlab

uiIsMatLab = false;
LIC = license('inuse');
for elem = 1:numel(LIC)
    envStr = LIC(elem).feature;
%    if strcmpi(envStr,'matlab')
    if ~isempty(strfind(lower(envStr),'matlab'))
        uiIsMatLab = true;
        break
    end
end
end % of function uiIsMatLab

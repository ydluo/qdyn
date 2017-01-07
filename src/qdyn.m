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
%		MESHDIM	dimension of the problem:
%			0 = spring-block system
%			1 = 1D fault in a 2D elastic medium
%			2 = 2D fault in a 3D medium
%		MU 	shear modulus (Pa)
%		LAM 	elastic modulus LAMBDA for 3D simulations (Pa)
%		VS 	shear wave velocity (m/s). If VS=0 radiation damping is turned off
%		L 	fault length if MESHDIM=1
%		    	stiffness is MU/L if MESHDIM=0
%		FINITE	boundary conditions if MESHDIM=1
%			0 = fault is infinitely long but slip is spatially periodic with period L,
%			    loaded by steady displacement at distance W from the fault
%			1 = fault is infinitely long but only a segment of length L has
%			    rate-and-state friction, the rest has steady slip. If you get the
%			    error message ???finite kernel is too small???, create a larger kernel file 
%			    using the function TabKernelFiniteFlt.m, update the file name in
%			    subroutine init_kernel_2D of src/fault_stress.f90, and recompile
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
%		TH_SS	reference steady-state state (default: TH_SS=DC/V_SS)
%		RNS_LAW	type of rate-and-state friction law:
%			0 = original
%			1 = with cutt-off velocities
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
%		N	number of fault elements if MESHDIM=1
%		NX	number of fault elements along-strike in 3D
%		NW 	number of fault elements along-dip in 3D
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
%		EXEC_PATH path to the Fortran qdyn executable
%       NPROCS = 1  % Default for serial runs
%              > 1  % MPI runs, change MPI_parallel=.true. in constants.f90
% OUTPUTS 	pars	structure containing the parameters listed above, and:
%			X,Y,Z = fault coordinates
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
% AUTHOR	Jean-Paul Ampuero	ampuero@gps.caltech.edu
% MODIFIED by Yingdi LUO        luoyd@gps.caltech.edu
% Last Mod 11/11/2014

function [pars,ot,ox] = qdyn(mode,varargin)

NPROCS = 1; % Default serial, no-MPI.

% NOTE ON VARIABLE NAMING CONVENTIONS
%	lower_case 	= local variables
%	UPPER_CASE 	= variables that will be wrapped into output structure 'pars'

%-- useful units
day = 60*60*24;
month = 30*day;
year = 365*day;

%--------- DEFAULT PARAMETERS ------------------------------------

MESHDIM=1;
NEQS=2;
NAME ='';	% title for the simulation

scriptName = mfilename('fullpath');
EXEC_PATH = fileparts(scriptName);	% default is same directory as qdyn.m

%-- medium
L= 2e3; 	% fault length (L scales the stiffness for the spring-block case)
FINITE=0;	% along strike bcs: 1=fixed-length asperity surrounded by steady creep, 0=periodic
W= 50e3;   	% out-of-plane dimension (ignored in spring-block)
MU = 30e9;	% shear modulus
LAM = 30e9;
VS = 3000; 	% shear wave velocity (if VS=0: turn off radiation damping)
V_TH= 1e-5; % threshold velocity for seismic events;

%-- numerical settings
N=1024/2; 	% number of grid cells
NX=100/2;
NW=10;
DW=1e3;
DIP_W=30.0;

Z_CORNER=-50e3;
IC=1;         %output ot coordinate
TMAX = 6*month;  % total simulation time
NSTOP = 0;	% stop at (0) tmax, (1) end of localization or (2) max slip rate
DTTRY = 1e2;   % first trial timestep
DTMAX = 0;	% maximum timestep (0=unrestricted)
ACC = 1e-7;     % solver accuracy
NXOUT = 8;	% space stride (cells) for snapshot outputs
NXOUT_DYN = 1;	% space stride (cells) for dynamic snapshot outputs
NTOUT = 100; 	% time stride (iterations) for snapshot outputs
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
TH_SS = DC/V_SS;
THETA_LAW = 1;
RNS_LAW = 0;
SIGMA_CPL = 0;
CO = 0;
%-- initial conditions
SIGMA=1e8;
V_0=V_SS ; 
TH_0=TH_SS;
V1=0.01;
V2=1e-7;
BRANCH='.false.';
STRIKE=0;
XC=0; % For branching fault, Junction location. 
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

pathstr = fileparts(mfilename('fullpath'));

% Override with inputs
Parse_Inputs(varargin{:});

X = [];
% fault mesh
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

TH_SS = DC./V_SS;

% wrap UPPER CASE variables in parameter structure fields with the same name 
fpars = who;
for k= find( strcmp(fpars,upper(fpars)) )' ,
  pars.(fpars{k}) = eval(fpars{k}) ;
end


switch mode

 case 'set', 
   ot=[];
   ox=[];
   return % pars = qdyn('set',...)  --> do not compute, exit here
  
 case 'read',
   pars = read_qdyn_in(NAME);
   pars.NAME = NAME;
   [pars.N,pars.FINITE] = read_qdyn_h(NAME);
   [ot,ox]= read_qdyn_out(NAME);

 case {'run', 'write'},

%     % recompile if qdyn.h must change
%     [n,finite] = read_qdyn_h(fullfile(pathstr,'qdyn'));
%     if (N~=n | FINITE~=finite )
%       write_qdyn_h(N,FINITE,pathstr);
%       cmd = ['cd ' pathstr '; make qdyn'];
%       system(cmd);
%     end
    
    % make vectors if constants
   DW(1:NW) =DW;  
   DIP_W(1:NW) =DIP_W;
   A(1:N)   =A;
   B(1:N)   =B;
   DC(1:N)  =DC;
   V_0(1:N) =V_0;
   TH_0(1:N)=TH_0;
   V1(1:N) =V1;
   V2(1:N) =V2;
   SIGMA(1:N) =SIGMA;
   MU_SS(1:N)=MU_SS;
   V_SS(1:N)=V_SS;
   IOT(1:N)=IOT;
   IASP(1:N)=IASP;
   CO(1:N)=CO;
    
    % For branching faults.
   if strcmp(BRANCH,'.true.')
     fid=fopen('qdyn_branch.in','w');    
     r=sqrt(X.^2+Y.^2);
     X=r.*cosd(STRIKE)+XC;
     Y=r.*sind(STRIKE)+YC;
     
     fprintf(fid,'%d\n',N);
     fprintf(fid,'%15.6f\n',DIP_W(1));
     fprintf(fid,'%20.6f %20.6f\n',LAM,MU);
     fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n',...
      [X(:),Y(:),Z(:),SIGMA(:),V_0(:),TH_0(:),A(:),B(:),DC(:),V1(:),V2(:),MU_SS(:),V_SS(:),CO(:)]');
     fclose(fid);
        ot = 0;
        ox = 0;
     return;
   end;
  %% Station
  if (MESHDIM==2);
   fids=fopen('stations.dat','w');
   nsta=1; % For now one station but later will be extended to more stations
   fprintf(fids,'%d\n',nsta);
   fprintf(fids,'%.15g %.15g %.15g\n',X(IC),Y(IC),Z(IC));
   fclose(fids);
  end;
    % export qdyn.in
  if (NPROCS>1); % MPI parallel 
   % Defining nwLocal 
   nwLocal(1:NPROCS)=floor(NW/NPROCS);
   % In case NW/NPRCOS is not integer. Leaving the rest the last processor
   nwLocal(NPROCS) = mod(NW,NPROCS) + nwLocal(NPROCS);
   nnLocal = 0;
   for iproc=0:NPROCS-1  
    iprocstr = num2str(sprintf('%06i',iproc));
    filename = ['qdyn' iprocstr '.in'];
    fid=fopen(filename,'w');
    fprintf(fid,'%u     meshdim\n' , MESHDIM); 
    if SIGMA_CPL == 1
      NEQS = 3;
    end
    if MESHDIM == 2;
      DWnnlocal=[];DIP_Wnnlocal=[];
      fprintf(1, 'MESHDIM = %d\n', MESHDIM); %JPA should "1" be "fid" instead?
      fprintf(fid,'%u %u     NX, NW\n' , NX, nwLocal(iproc+1));      
      fprintf(fid,'%.15g %.15g  %.15g      L, W, Z_CORNER\n', L, W, Z_CORNER);
      DWnnlocal = DW(nnLocal/NX+1:nwLocal(iproc+1)+nnLocal/NX);
      DIP_Wnnlocal = DIP_W(nnLocal/NX+1:nwLocal(iproc+1)+nnLocal/NX);
      fprintf(fid,'%.15g %.15g \n', [DWnnlocal(:),DIP_Wnnlocal(:)]');
    else  
      fprintf(fid,'%u     NN\n' , N);      
      fprintf(fid,'%.15g %.15g      L, W\n', L, W);
    end
    
    if MESHDIM == 1;
        fprintf(fid,'%u   finite\n', FINITE);
    end   
    
    fprintf(fid,'%u   itheta_law\n', THETA_LAW);
    fprintf(fid,'%u   i_rns_law\n', RNS_LAW);
    fprintf(fid,'%u   i_sigma_cpl\n', SIGMA_CPL);    
    fprintf(fid,'%u   n_equations\n', NEQS);
    fprintf(fid,'%u %u %u %u %u %u  ntout, nt_coord, nxout, nxout_DYN, ox_SEQ, ox_DYN\n', NTOUT,IC,NXOUT,NXOUT_DYN,OX_SEQ,OX_DYN);     
    fprintf(fid,'%.15g %.15g %.15g %.15g    beta, smu, lambda, v_th\n', VS, MU, LAM, V_TH);
    fprintf(fid,'%.15g %.15g    Tper, Aper\n',TPER,APER);
    fprintf(fid,'%.15g %.15g %.15g %.15g    dt_try, dtmax, tmax, accuracy\n',DTTRY,DTMAX,TMAX,ACC);
    fprintf(fid,'%u   nstop\n',NSTOP);
    fprintf(fid,'%u %u  DYN_FLAG, DYN_SKIP\n',DYN_FLAG,DYN_SKIP);
    fprintf(fid,'%.15g %.15g %.15g    M0, DYN_th_on, DYN_th_off\n', DYN_M,DYN_TH_ON,DYN_TH_OFF);

    for wloc=1:nwLocal(iproc+1)
     for xloc=1:NX
      iloc = xloc + (wloc-1)*NX + nnLocal;
      fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n',...
        SIGMA(iloc),V_0(iloc),TH_0(iloc),A(iloc),B(iloc),DC(iloc),V1(iloc),V2(iloc),MU_SS(iloc),V_SS(iloc),IOT(iloc),IASP(iloc),CO(iloc));
     end;
    end;
  
    for wloc=1:nwLocal(iproc+1)
     for xloc=1:NX
      iloc = xloc + (wloc-1)*NX + nnLocal;
      fprintf(fid,'%.15g %.15g %.15g %.15g\n',...
          X(iloc),Y(iloc),Z(iloc),DIP(iloc));
     end;
    end;
    
    % next processor
    nnLocal=NX*nwLocal(iproc+1) + nnLocal;
    % hold on
    %scatter(X(1+nnLocal:iloc),Z(1+nnLocal:iloc),[],Z(1+nnLocal:iloc))
    fclose(fid);    
   end;
  end;

  fid=fopen('qdyn.in','w');
   fprintf(fid,'%u     meshdim\n' , MESHDIM); 
   if SIGMA_CPL == 1
       NEQS = 3;
   end
   if MESHDIM == 2;
       fprintf(1, 'MESHDIM = %d\n', MESHDIM); %JPA should "1" be "fid" instead?
       fprintf(fid,'%u %u     NX, NW\n' , NX, NW);      
       fprintf(fid,'%.15g %.15g  %.15g      L, W, Z_CORNER\n', L, W, Z_CORNER);
       fprintf(fid,'%.15g %.15g \n', [DW(:), DIP_W(:)]');
   else  
       fprintf(fid,'%u     NN\n' , N);      
       fprintf(fid,'%.15g %.15g      L, W\n', L, W);
   end
    
   if MESHDIM == 1;
         fprintf(fid,'%u   finite\n', FINITE);
   end   
    
   fprintf(fid,'%u   itheta_law\n', THETA_LAW);
   fprintf(fid,'%u   i_rns_law\n', RNS_LAW);
   fprintf(fid,'%u   i_sigma_cpl\n', SIGMA_CPL);    
   fprintf(fid,'%u   n_equations\n', NEQS);
   fprintf(fid,'%u %u %u %u %u %u  ntout, nt_coord, nxout, nxout_DYN, ox_SEQ, ox_DYN\n', NTOUT,IC,NXOUT,NXOUT_DYN,OX_SEQ,OX_DYN);     
   fprintf(fid,'%.15g %.15g %.15g %.15g    beta, smu, lambda, v_th\n', VS, MU, LAM, V_TH);
   fprintf(fid,'%.15g %.15g    Tper, Aper\n',TPER,APER);
   fprintf(fid,'%.15g %.15g %.15g %.15g    dt_try, dtmax, tmax, accuracy\n',DTTRY,DTMAX,TMAX,ACC);
   fprintf(fid,'%u   nstop\n',NSTOP);
   fprintf(fid,'%u %u  DYN_FLAG, DYN_SKIP\n',DYN_FLAG,DYN_SKIP);
   fprintf(fid,'%.15g %.15g %.15g    M0, DYN_th_on, DYN_th_off\n', DYN_M,DYN_TH_ON,DYN_TH_OFF);

          
   fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n',...
      [SIGMA(:),V_0(:),TH_0(:),A(:),B(:),DC(:),V1(:),V2(:),MU_SS(:),V_SS(:),IOT(:),IASP(:),CO(:)]');
 
  fclose(fid);
    
    if strcmp(mode, 'write')
        ot = 0;
        ox = 0;
        return;
    end
    
%    Solve
    if (NPROCS==1) 
       status = system([EXEC_PATH filesep 'qdyn']);
    else
       status = system(['mpirun -np ' num2str(NPROCS) ' ' EXEC_PATH filesep 'qdyn']);
    end
%   rename input and output files
    if length(NAME)
      movefile('fort.18',[NAME '.ot']); 
      movefile('fort.19',[NAME '.ox']); 
      copyfile('qdyn.in',[NAME '.in']); 
      copyfile(fullfile(pathstr,'qdyn.h') ,[NAME '.h']); 
    end
    % output
    if (NPROCS>1); % MPI parallel
      [ot,ox]= read_qdyn_out_mpi(NAME);
    else
      [ot,ox]= read_qdyn_out(NAME);
    end

  otherwise,
    error('mode must be: set, read or run')

end

%-----------
% PARSE_INPUTS sets variables in the caller function according to inputs.
%	By convention the variables to set have UPPER CASE names. 
%
% WARNING: no checks yet 
% Field names in parsin and Property's should match field names in pars
% or else the default value will be taken.
function Parse_Inputs(varargin)

if isempty(varargin), return, end

% 1. Override defaults by input structure
if isstruct(varargin{1}) 
 % process input parameter structure
  parsin = varargin{1};
  fparsin = fieldnames(parsin);
  for k=1:length(fparsin),
    assignin('caller', upper(fparsin{k}), parsin.(fparsin{k}) );
  end
 % separate the input Property/Value pairs
  varargin(1)=[];
end

% 2. Override defaults and previously set parameters by Property/Value pairs
for k=2:2:length(varargin),
  assignin('caller', upper(varargin{k-1}), varargin{k} );
end


%-----------
% function [N,FINITE] = read_qdyn_h(name)
% if ~exist('name','var') || ~length(name), name = 'qdyn'; end
% name = [name '.h'];
% fid=fopen(name);
% 
% rline=fgetl(fid);
% inn=strfind(rline,'nn=');
% N=sscanf(rline(inn:end),'nn=%u');
% 
% rline=fgetl(fid);
% inn=strfind(rline,'finite=');
% FINITE=sscanf(rline(inn:end),'finite=%u');
% 
% fclose(fid);

%-----------
% function write_qdyn_h(N,FINITE,pathstr)
% fid=fopen(fullfile(pathstr,'qdyn.h'),'w');
% fprintf(fid,'      parameter(nn=%u)\n',N);
% fprintf(fid,'      parameter(finite=%u)\n',FINITE);
% nnfft=(FINITE+1)*N;
% fprintf(fid,'      parameter(nnfft=%u)\n',nnfft);
% fprintf(fid,'      parameter(nwfft=%u)\n',2+ceil(sqrt(nnfft/2)) );
% fclose(fid);

%-----------
function pars = read_qdyn_in(name)

if ~exist('name','var') || ~length(name), name = 'qdyn'; end
name = [name '.in'];

fid=fopen(name);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
MU_SS = rdat(1);
V_SS  = rdat(2);
THETA_LAW = rdat(3);
SIGMA = rdat(4);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
MU = rdat(1);
VS = rdat(2);
W  = rdat(3);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
L     = rdat(1);
NSTOP = rdat(2);
TMAX  = rdat(3);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
TPER  = rdat(1);
APER  = rdat(2);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
NXOUT = rdat(1);
NTOUT = rdat(2);
rline=fgetl(fid); rdat = sscanf(rline,'%f');
DTTRY = rdat(1);
DTMAX = rdat(2);
ACC   = rdat(3);
fclose(fid);
[X,A,B,DC,V1,V2,V_0,TH_0,]=textread(name,'','headerlines',6);

% wrap UPPER CASE variables in parameter structure fields with the same name 
fpars = who;
for k= find( strcmp(fpars,upper(fpars)) )' ,
  pars.(fpars{k}) = eval(fpars{k}) ;
end


% read outputs from qdyn.f
function [ot,ox] = read_qdyn_out_mpi(name)

if exist('name','var') && length(name)
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
  ox.vd= cosa(:,:,7);
  ox.dtau = cosa(:,:,8);
  ox.dtaud = cosa(:,:,9);
  ox.d = cosa(:,:,10);
  ox.sigma = cosa(:,:,11);

%-----------
% read outputs from qdyn.f
function [ot,ox] = read_qdyn_out(name)

if exist('name','var') && length(name)
  namet = [name '.ot'];
  namex = [name '.ox'];
else
  namet = 'fort.18';
  namex = 'fort.19';
end

if uimatlab
  % time series
  [ot.t, ot.locl, ot.cl, ot.p, ot.pdot, ...
   ot.vc, ot.thc, ot.omc, ot.tauc, ot.dc, ...
   ot.xm, ot.v, ot.th, ot.om, ot.tau, ot.d, ot.sigma ] = ...
    textread(namet,'','headerlines',6);
  
  % snapshots
  fid=fopen(namex);
  NSX=fscanf(fid,'# nx=%u');
  fclose(fid);
  cosa = textread(namex,'','commentstyle','shell');
  ncosa = size(cosa);
  NST=ncosa(1)/NSX;
  cosa=reshape(cosa,NSX,NST,ncosa(2));
  ox.x = cosa(:,1,1);
  ox.t = cosa(1,:,2)';
  ox.v = cosa(:,:,3);
  ox.th= cosa(:,:,4);
  ox.vd= cosa(:,:,5);
  ox.dtau = cosa(:,:,6);
  ox.dtaud = cosa(:,:,7);
  ox.d = cosa(:,:,8);
  ox.sigma = cosa(:,:,9);

else 
  [ot,ox] = read_qdyn_out_Octave(namet,namex)
end

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

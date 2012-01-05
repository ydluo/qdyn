
% QDYN		Quasi-dynamic earthquake cycles on a 1D fault embedded in a 2D elastic medium
%               This is a Matlab wrapper for the Fortran code qdyn.f
%               Friction Law with cut-off Velocities by Okubo, velocity
%               weakening at low slip_v and strengthening at high slip_v
%		Features:
% 		+ rate-and-state friction 
%   			Mu = Muss + a*ln(V/Vss) + b*ln(Theta/Thetass)
%	  	  with ageing law:
%   			dTheta/dT = 1-V*Theta/Dc
%		  or slip law
%   			dTheta/dT = -V*Theta/Dc * log(V*Theta/Dc)
%                 or ageing law in the self-accelerating approximation:
%   			dTheta/dT = -V*Theta/Dc
%		+ spatially non-uniform a,b,Dc,v(0),theta(0)
%		+ quasistatic stress balance with radiation damping, no inertia, no elastodynamics
%		+ two possible boundary conditions: 
%			. the fault is periodic along-strike and is loaded by steady displacement
%			  at a fault-normal distance W
%			  (crustal plane model, mimics a finite seismogenic depth W)
%		  	. the fault area governed by rate-and-state friction has a 
%			  finite length and is loaded by steady sliding on the rest of the fault
%
% SYNTAX	[pars,ot,ox] = qdyn(mode,[parsin],['Property',Value,...])
%
% INPUTS 	mode	'set'	gives the default parameter structure (pars), 
%				overrides by fields present in structure parsin
%				or by Property/Value pairs
%			'run'	sets parameters and runs a simulation
%			'read' 	reads parameters and outputs from a previous simulation
%		parsin	parameter structure to override the default parameters 
%		'Prop' 	property to be set, a fieldname of the parameter structure
%		Value   of the above property, overrides default and parsin
%		
%		These are the parameters that can be set through 'parsin' or 'Prop/Value' pairs:
%
%		L = fault length (L scales the stiffness for the spring-block case)
%		FINITE = boundary conditions: 
%			0 = periodic along-strike, steady loading at distance W from the fault line
%			1 = rate-and-state fault segment of finite length (L) surrounded by steady slip
%		W  = out-of-plane dimension (ignored if FINITE=1)
%		MU = shear modulus
%		VS = shear wave velocity. If VS=0 radiation damping is turned off
%		N  = number of fault nodes (elements)
%		TMAX = total simulation time (in seconds)
%		NSTOP = stop at (0) t=TMAX, (1) end of localization or (2) first slip rate peak
%		DTTRY = first trial timestep (in seconds)
%		DTMAX = maximum timestep (0=unrestricted)
%		ACC = solver accuracy
%		NXOUT = spatial interval (number of nodes) for snapshot outputs
%		NTOUT = temporal interval (number of iterations) for snapshot outputs
%		A  = amplitude of direct effect in rate-and-state friction 
%		B  = amplitude of evolution effect in rate-and-state friction
%		DC = characteristic slip in rate-and-state friction
%		MU_SS = reference steady-state friction coefficient
%		V_SS = reference steady-state slip velocity
%		TH_SS = reference steady-state state (normally TH_SS=DC/V_SS)
%		THETA_LAW = evolution law for the state variable:
%			0 = ageing in the "no-healing" approximation
%			1 = ageing law
%			2 = slip law
%		SIGMA = effective normal stress
%		V_0 = initial slip velocity
%		TH_0 = initial state 
%		APER = amplitude of additional periodic loading (in Pa)
%		TPER = period of additional periodic loading (in s)
%
% OUTPUTS 	pars	structure containing parameters, see documentation of qdyn.f
%		ot	structure containing time series outputs 
%			at the point of maximum slip rate
%			ot.t	output times
%			ot.locl	localization length (distance between stressing rate maxima)
%			ot.cl	crack length (distance between slip rate maxima)
%			ot.p	seismic potency
%			ot.pdot	seismic potency rate
%			ot.xm	location of maximum slip rate 
%			ot.v	maximum slip rate 
%			ot.th	state variable theta at xm
%			ot.om 	slip_rate*theta/dc at xm
%			ot.tau	stress at xm
%			ot.d	slip at xm
%			ot.vc	slip rate at center
%			ot.thc	state variable theta at center
%			ot.omc 	slip_rate*theta/dc at center
%			ot.tauc	stress at center
%			ot.dc	slip at center
%		ox	structure containing snapshot outputs (x,t)
%			ox.x	fault coordinates 
%			ox.t 	output times
%			ox.v	slip rate (x,t)
%			ox.th	state variable theta(x,t)
%			ox.vd	slip acceleration
%			ox.dtau stress (-initial)
%			ox.dtaud stress rate
%			ox.d 	slip
%
% EXAMPLE	A run with initial velocity slightly above the default steady state:
%			p = qdyn('set');
%			[p,ot,ox] = qdyn('run','V_0',1.01*p.V_SS);
%			semilogy(ot.t,ot.v)
%
% SEE ALSO	qdyn_example1.m
%
% AUTHOR	Jean-Paul Ampuero	ampuero@gps.caltech.edu
% LAST MODIF	May 2010

function [pars,ot,ox] = qdyn(mode,varargin)


% NOTE ON VARIABLE NAMING CONVENTIONS
%	lower_case 	= local variables
%	UPPER_CASE 	= variables that will be wrapped into output structure 'pars'

%--------- DEFAULT PARAMETERS ------------------------------------

MESHKIND=0;
KERNELKIND=0;
NEQS=2;


%-- useful units
day = 60*60*24;
month = 30*day;
year = 365*day;

NAME ='';	% title for the simulation

%-- medium
L= 2e3; 	% fault length (L scales the stiffness for the spring-block case)
FINITE=0;	% along strike bcs: 1=fixed-length asperity surrounded by steady creep, 0=periodic
W= 50e3;   	% out-of-plane dimension (ignored in spring-block)
MU = 30e9;	% shear modulus
VS = 3000; 	% shear wave velocity (if VS=0: turn off radiation damping)

%-- numerical settings
N=1024; 	% number of grid cells
IC=512;         %output ot coordinate
TMAX = 6*month;  % total simulation time
NSTOP = 0;	% stop at (0) tmax, (1) end of localization or (2) max slip rate
DTTRY = 1e2;   % first trial timestep
DTMAX = 0;	% maximum timestep (0=unrestricted)
ACC = 1e-7;     % solver accuracy
NXOUT = 8;	% space stride (cells) for snapshot outputs
NTOUT = 100; 	% time stride (iterations) for snapshot outputs

%-- friction
A = 0.9e-2; 
B = 1e-2; 
DC = 4e-4;	
MU_SS = 0.6;
V_SS = 1e-9;
TH_SS = DC/V_SS;
THETA_LAW = 1;

%-- initial conditions
SIGMA=1e8;
V_0=V_SS ; 
TH_0=TH_SS; 
V1=0.01;
V2=1e-7;


%-- periodic loading 
APER = 0;
TPER = 1*year; 

%--------- INITIALIZE -------------------------------------

pathstr = fileparts(mfilename('fullpath'));

% Override with inputs
Parse_Inputs(varargin{:});

X = (-N/2+0.5:N/2-0.5) *L/N; % fault coordinates
TH_SS = DC/V_SS;

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

  case 'run',

%     % recompile if qdyn.h must change
%     [n,finite] = read_qdyn_h(fullfile(pathstr,'qdyn'));
%     if (N~=n | FINITE~=finite )
%       write_qdyn_h(N,FINITE,pathstr);
%       cmd = ['cd ' pathstr '; make qdyn'];
%       system(cmd);
%     end
    
    % make vectors if constants
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
 
    % export qdyn.in
    fid=fopen('qdyn.in','w');
    if MESHKIND == 0;
      fprintf(fid,'%u     meshkind\n' , MESHKIND);
      fprintf(fid,'%u     NN\n' , N);      
      fprintf(fid,'%.15g %.15g      L, W\n', L, W);
      if KERNELKIND ==0;
          fprintf(fid,'%u   kernelkind\n', KERNELKIND);
          fprintf(fid,'%u   finite\n', FINITE);
          fprintf(fid,'%u   itheta_law\n', THETA_LAW);
          fprintf(fid,'%u   n_equations\n', NEQS);
          fprintf(fid,'%u %u %u   ntout, nt_coord, nxout\n', NTOUT,IC,NXOUT);     
          fprintf(fid,'%.15g %.15g   beta, smu\n', VS, MU);
          fprintf(fid,'%.15g %.15g    Tper, Aper\n',TPER,APER);
          fprintf(fid,'%.15g %.15g %.15g %.15g    dt_try, dtmax, tmax, accuracy\n',DTTRY,DTMAX,TMAX,ACC);
          fprintf(fid,'%u   nstop\n',NSTOP);
          
          fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n',[SIGMA(:),V_0(:),TH_0(:),A(:),B(:),DC(:),V1(:),V2(:),MU_SS(:),V_SS(:)]');
    
      end
    end
    fclose(fid);
    
    % solve
%     status = system('~/qdyn_svn/trunk/src/qdyn');
    status = system('~/bin/qdyn');
    % rename input and output files
    if length(NAME)
      movefile('fort.18',[NAME '.ot']); 
      movefile('fort.19',[NAME '.ox']); 
      copyfile('qdyn.in',[NAME '.in']); 
      copyfile(fullfile(pathstr,'qdyn.h') ,[NAME '.h']); 
    end

    % output
    [ot,ox]= read_qdyn_out(NAME);

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
[X,A,B,DC,V1,V2,V_0,TH_0]=textread(name,'','headerlines',6);

% wrap UPPER CASE variables in parameter structure fields with the same name 
fpars = who;
for k= find( strcmp(fpars,upper(fpars)) )' ,
  pars.(fpars{k}) = eval(fpars{k}) ;
end

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

% time series
[ot.t, ot.locl, ot.cl, ot.p, ot.pdot, ...
 ot.vc, ot.thc, ot.omc, ot.tauc, ot.dc, ...
 ot.xm, ot.v, ot.th, ot.om, ot.tau, ot.d ] = ...
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


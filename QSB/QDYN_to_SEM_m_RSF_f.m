clear;
clc;

t_dyn = 300;        %time elipsed in SPECFEM simu

x_off=-256e3;
y_off=0;
z_off=0;

%pars
for istart=1:1:1
mname='nodesonfault'; %SEM mesh
% wtname='input_file.txt';  %SEM input file
NGLL=4;
itp=1;   %interpfactor


rn=200;
rsn=3;
vmax=1.;
%c_depth=3e3; % Cohesion cutoff depth in m
%coh=4e6;    %cohension in Pa/ default: 4e6
%forced rupture zone
hxr=0.5;  %hypocenter relative location along-strike (0,1) left to right
hzr=10/50; %hypocenter relative depth (0,1) surface to bottom
rrpt=10e3; %radius of forced rupture zone in m
vrr=2e3;  %forced rupture prop speed in m/s
hypo_auto = 1;  %auto detect hypocenter location

idummy=1;
fdummy=1.;



%load mesh
disp(['Loading SEM mesh :  ',mname,'...']);
xyz=load(mname);
q.X=xyz(:,3);
q.Y=xyz(:,4);
q.Z=xyz(:,5);
q.N=numel(xyz(:,3));
%find along strike length
q.NZ=0;
while q.X(q.NZ+1) == q.X(1)
    q.NZ=q.NZ+1;
end
q.NX=q.N/q.NZ;


disp('Done');

disp('Loading qdyn.in ...');
p=Qdyn_read_in();
disp('Done');


%load QDYN output
rname=['fort.' num2str(19998+istart*3)];
disp(['Loading QDYN output :  ',rname,'...']);
o = Qdyn_read_ox_seq(rname);
disp('Done');

p.X=o.X+x_off;
p.Y=o.Y+y_off;
p.Z=o.Z+z_off;


%Reorganize coordinate
%locate hypocenter
hx=p.L*hxr;
hy=mean(p.Y);
hz=p.Z_CORNER*hzr;

if hypo_auto == 1
    [vmax_i,loc_i]=max(o.v);
    hx=p.X(loc_i);
    hy=p.Y(loc_i);
    hz=p.Z(loc_i);
end
    

co=zeros(q.N,1);
tau=zeros(q.N,1);

disp('Matching Values ...');


for i=1:q.N
    
    [cc, id] = min((p.X(:)-q.X(i)).^2 + (p.Z(:)-q.Z(i)).^2);
    
    q.SIGMA(i)=p.SIGMA(id);
    q.A(i)=p.A(id);
    q.B(i)=p.B(id);
    q.DC(i)=p.DC(id);
    q.V0(i)=p.V_SS(id);
    q.f0(i)=p.MU_SS(id);
    q.V_ini(i)=o.v(id);
    q.TH_ini(i)=o.th(id);
    co(i)=p.CO(id);
    tau(i)=q.SIGMA(i)*(q.f0(i)+q.A(i)*log(q.V_ini(i)/q.V0(i))+q.B(i)*log(q.TH_ini(i)*q.V0(i)/q.DC(i)));
    
%    if abs(q.Z(i)) <= c_depth
%        co(i)=coh;
%    end

    if mod(i,ceil(q.N/100)) == 0
        disp([num2str(i/ceil(q.N/100)) '% Complete'])
    end
end

disp('Done');




fid=fopen('rsf_hete_input_file.txt','w');

fprintf(fid,'%u %u %E %E \n' ,q.NX, q.NZ, max(q.X)-min(q.X), max(abs(q.Z))-min(abs(q.Z))); 
disp('Start write SEM input:...');
for id=1:q.N
    fprintf(fid,'%E %E %E %E %E %E %E %E %E %E %E %E %E\n' , ...
      q.X(id),q.Z(id),...
      q.SIGMA(id),tau(id),0.,q.V0(id),q.f0(id),...
      q.A(id),q.B(id),q.DC(id),q.V_ini(id),q.TH_ini(id),co(id));
%    if mod(id,100) == 1
%          disp(['Wrote ',num2str(id),' of ',num2str(q.N)])
%    end

    if mod(id,ceil(q.N/100)) == 0
        disp([num2str(id/ceil(q.N/100)) '% Complete'])
    end
end
disp(['Generated SEM input :  input_file.txt']);


fclose(fid);

if ~(exist('timestamp.txt'))
    fid=fopen('timestamp.txt','w');
    fprintf(fid,'%E %E',0,o.t);
    fclose(fid);
    display(['t0 = ' num2str(0) '     t1 = ' num2str(o.t)])
else
    fid=fopen('timestamp.txt','r');
    rline=fgetl(fid); rdat = sscanf(rline,'%E');
    t0=rdat(1);
    t1=rdat(2);
    display(['Previous t0 = ' num2str(t0) '     t1 = ' num2str(t1)])
    fclose(fid);
    fid=fopen('timestamp.txt','w');
    fprintf(fid,'%E %E',t1+t_dyn,t1+t_dyn+o.t);
    display(['t0 = ' num2str(t1+t_dyn) '     t1 = ' num2str(t1+t_dyn+o.t)])
    fclose(fid);
    
end



end

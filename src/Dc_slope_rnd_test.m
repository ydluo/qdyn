% Test of Dc slope controlling rupture

clear;
clc;

%------------------------------
rand(1,floor(sum(100*clock)));
%------------------------------


L = 128e3;
N = 1024;
W = 15e3;

FINITE = 1;

twm2 = 4000;     %simu time in years 
      % slope left dDc/dx normalized by critical DDCC = sigma(b-a)/mu

bin_DC_mean = [0.01:0.005:0.1 0.012:0.005:0.1];		%mean DC
bin_dcsigma = [0.25 0.5 0.125];		%Dv
bin_col_l = [125:125:1000 1250:250:5000 6e3:1e3:10e3 12e3:2e3:30e3];	%r of col 

%bin_DC_mean = [0.025];		%mean DC
%bin_dcsigma = [0.25];		%Dv
%bin_col_l = [5000];	%r of col 



DC_min = 0.01;		%min_DC




year = 3600*24*365;      
      
      
%-----------


tardir = ' luoyd@pacha:/export/raid1/luoyd/qdyn_results_DC_test_rnd/';
filename_sum = 'Dc_test_rnd_finiti.txt';
event_type = {'Undefined','Steady-Slip','Irreg. Aseismic','Chara. Aseismic','Irreg. Seismic','Chara. Seismic'};

t_pk_dist = 1800;        % min peak distance in seconds

v_th = 0.1; %threshold for dynamic events;
v_th_ss = 1.01; 
t_rec_th = 0.2; %threshold for irregular events; t_rec_min < t_rec_max * (1-t_rec_th)
num_cc = 5;     % number of bins used to decide if a event is cont. dec.
t_end_p = 0.5;  %portion of twm used for analysiscd 


ii_count_glb = 0;


% dia.HC = zeros(ttl_sim,1);
% dia.L_slp = zeros(ttl_sim,1);
% dia.R_slp = zeros(ttl_sim,1);
% dia.L_rup = zeros(ttl_sim,1);
% dia.R_rup = zeros(ttl_sim,1);
% dia.Dc_rup = zeros(ttl_sim,1);
% dia.Lc_rup = zeros(ttl_sim,1);
% dia.vmax = zeros(ttl_sim,1);
% dia.event_type = zeros(ttl_sim,1);


for iibin_dcsigma = 1:1:numel(bin_dcsigma)
    dcsigma = bin_dcsigma(iibin_dcsigma);
    
    for iiDC_mean = 1:1:numel(bin_DC_mean)
        DC_mean = bin_DC_mean(iiDC_mean);
        
        for iicol_l = 1:1:numel(bin_col_l)
            col_l = bin_col_l(iicol_l);
                    
           
            ii_count_glb = ii_count_glb + 1;
            
            DC_v = (exp((dcsigma)^2)-1)*DC_mean^2;
            dcmu = log(DC_mean^2/sqrt(DC_v+DC_mean^2));
            
            filename = ['Dc_test_rnd_finite_dcsigma' num2str(dcsigma) '_DCmean' num2str(DC_mean) ...
                '_coll' num2str(col_l) '.mat'];
            
            disp(filename);

            year = 3600*24*365;
            p = qdyn('set');
            p.NSTOP = 0;        %stop at v_th
            p.TMAX = 100000*year;       % stop at v = v_th = tmax


            co = 4e6;  %cohesion
            co_limit = 3e3;  %first X m to appy cohesion 

            p.RNS_LAW=0;
            p.MESHDIM=2;      %FFT enabled
            p.THETA_LAW=1;
            p.MU=40e9;
            p.LAM=40e9;
            p.MU_SS=0.6;
            %p.SIGMA=100e6;
            p.V_SS=0.01/year;
            p.V2=100.;              %no cut off velocity
            p.V1=p.V2;

            p.OX_SEQ=1;
            p.OX_DYN=1;
            p.DYN_TH_ON = 0.1;
            p.DYN_TH_OFF = 0.1;
            %p.DC=0.3;



            dip0=90.;
            dw0=0.5e3/4;
            db=20e3;     %bottom of simulation zone (depth in m)
            aa0=0.01;    %p.A
            sigma0=75e6;        %sigma max
            nxout=1;       %snapshot output grid interval
            p.NXOUT_DYN=1;   %dynamic snapshot output grid interval
            p.W=db/sin(dip0/180.*pi);
            p.L=512e3;
            p.NW=ceil(p.W/dw0);
            p.NX=1024*4;
            p.VS=3000.;

            r_filter = ceil(col_l/dw0);

            p.N=p.NX*p.NW;
            p.DW(1:p.NW)=p.W/p.NW;     %deep to shallow
            p.DIP_W(1:p.NW)=dip0;      %deep to shallow 
            %p.Z_CORNER=-db+.5*p.DW(1)*sin(p.DIP_W(1)/180.*pi);   %p.Z_CORNER at center of left-bottom cell

            p.Z_CORNER=-db;

            p.IC = p.N/2;

            p.A(1:p.NW)=aa0;
            dz0=dw0*sin(dip0/180.*pi);
            ba0=1.5;     %b/a at seismogenic zone
            %abmax=5;     %a/b at d3

            p.B(1:p.NW)=p.A(1:p.NW).*ba0;

            p.DC(1:p.NW)=0.3;

            p.SIGMA(1:p.NW)=sigma0;

            p.CO(1:ceil((db-co_limit)/dz0))=0;
            p.CO(ceil((db-co_limit)/dz0)+1:p.NW)=co;

            p.N=p.NX*p.NW;
            tmp_A=p.A;
            tmp_B=p.B;
            tmp_SIGMA=p.SIGMA;
            tmp_DC=p.DC;
            tmp_CO=p.CO;


            for i=1:p.NW
                p.X((i-1)*p.NX+1:i*p.NX) = linspace(0,p.L,p.NX);
                p.A((i-1)*p.NX+1:i*p.NX) = tmp_A(i);
                p.B((i-1)*p.NX+1:i*p.NX) = tmp_B(i);
                p.SIGMA((i-1)*p.NX+1:i*p.NX) = tmp_SIGMA(i);
                p.DC((i-1)*p.NX+1:i*p.NX) = tmp_DC(i);
                p.CO((i-1)*p.NX+1:i*p.NX) = tmp_CO(i);

            end


            for i=1:1:p.N
                p.IOT(i) = 0;
            end

            T = normrnd(log(DC_mean),dcsigma,p.NX+r_filter*2,p.NW+r_filter*2);
            h = fspecial('disk',r_filter);
            Tm = imfilter(T,h);
            TTm = Tm(r_filter+1:1:end-r_filter,r_filter+1:1:end-r_filter);
            pd = fitdist(reshape(TTm,[],1),'normal');
            TTm = pd.mu+(TTm-pd.mu)*dcsigma/pd.sigma;
            TTm = exp(TTm);
            p.DC = reshape(TTm,size(p.X));
            p.DC = max(p.DC,DC_min);


            twm=100000;         %warmup time in years
            ts=1000;    %simulation time in years
            p.ACC = 1e-10;
            Vdyn=2*mean(p.A.*p.SIGMA./p.MU.*p.VS);
            p.DYN_TH_ON=0.01;
            p.DYN_TH_OFF=0.01*0.5;

            %------------------------------
            Lb = min(p.MU.*p.DC./p.SIGMA./p.B)
            %Lnuc = 1.3774*Lb;
            %------------------------------


            p.IC=ceil(p.N/2);
            dw=p.W/p.NW;
            Lb_over_dw = Lb/dw;
            dx=p.L/p.NX;
            Lb_over_dx = Lb/dx;

            p = qdyn('set',p);

            p.V_0 = 1.01*p.V_SS ;


            p.NTOUT=1000;
            p.NXOUT=nxout;

            % p.DYN_FLAG=1;
            % p.DYN_M=10.^19;
            % p.DYN_SKIP = 1;
            
            
            p2 = qdyn('set');
            p2.B = p.B(p.N/2+1);
            p2.A = p.A(p.N/2+1);
            p2.L = L;
            p2.FINITE=FINITE;
            p2.N = N;
            p2.DC = p.DC(p.N/2+1:1:p.N/2+p2.N);
            p2.W = W;
            p2.SIGMA = p.SIGMA(p.N/2+1);
            p2.MU = p.MU;
            Lb_min = p2.MU*min(p2.DC)/p2.SIGMA/p2.B;
            Lb_mean = p2.MU*mean(p2.DC)/p2.SIGMA/p2.B;




            dx=p2.L/p2.N;
            Lb_min_over_dx = Lb_min/dx
            Lb_mean_over_dx = Lb_mean/dx

            p2.ACC = 1e-10;

            p2 = qdyn('set',p2);


            p2.TMAX = twm2 * year; % 
            p2.NTOUT=100;
            p2.NXOUT=1;
            p2.NSTOP=0;
            

            
            p2.V_0 = ( 1.+0.01*exp( -(p2.X/(p2.L*0.1)).^6 ) )*p2.V_SS ;
            p2.V_0 = p2.V_0/mean(p2.V_0)*p2.V_SS;

            [p2,ot,ox] = qdyn('run',p2) ;
            
            Vdyn=2*mean(p2.A.*p2.SIGMA./p2.MU.*p2.VS);

            filename_eps = [filename '.eps'];
            filename_eps_z = [filename '_zoom.eps'];

            
            
                  %------------ try to find event type

            ii_end = find(ot.t >= twm2*year*(1.0 - t_end_p) , 1);
            vend = ot.v(ii_end:1:numel(ot.v));
            tend = ot.t(ii_end:1:numel(ot.v));

            vend_max = max(vend);

            i_event_type = 1;

            isp = 0;
            if vend_max >= v_th * Vdyn 
                i_event_type = 6;   %dynamic event
                [pks,locs] = findpeaks(vend,tend,'MinPeakDistance',t_pk_dist,'MinPeakHeight',v_th * Vdyn * 0.1);
                t_rec = diff(locs);
                if min(t_rec) < max(t_rec) * (1.0 - t_rec_th)   %irregular seismic
                    i_event_type = 5;
                end
                if numel(t_rec)<=1          % only one or less trec is find is find
                    [pks,locs] = findpeaks(vend,tend,'MinPeakDistance',t_pk_dist,'MinPeakHeight',v_th * Vdyn * 0.001);
                    t_rec = diff(locs);
                    if numel(t_rec) > 1
                        i_event_type = 5;
                    else
                        isp = 1;
                    end

                end

            end

            if vend_max <= p2.V_SS*v_th_ss
                    i_event_type = 2;   %steady slip
            end

            if ((vend_max > p2.V_SS*v_th_ss) && (vend_max < v_th * Vdyn)) || (isp == 1)
                    i_event_type = 4;   % chara aseis
                [pks,locs] = findpeaks(vend,tend,'MinPeakDistance',t_pk_dist,'MinPeakHeight',vend_max * 0.01);
                t_rec = diff(locs);
                if min(t_rec) < max(t_rec) * (1.0 - t_rec_th)   %irregular aseismic
                    i_event_type = 3;
                end

                if (vend_max <= p2.V_SS*100)
                    ii_end_all = zeros(num_cc,1);
                    vend_max_all = zeros(num_cc-1,1);

                    for ii_cc = 1:1:num_cc
                        ii_end_all(ii_cc) = find(ot.t >= twm*year*(0.99999 - 0.5/num_cc*(num_cc-ii_cc)) , 1);
                    end

                    for ii_cc = 1:1:num_cc-1
                        vend_max_all(ii_cc) = max(ot.v(ii_end_all(ii_cc):1:ii_end_all(ii_cc+1)));
                    end

                    vend_max_all = fliplr(vend_max_all);

                    if issorted(vend_max_all)
                        i_event_type = 2;
                    end  

                    if (vend_max <= p2.V_SS*2)
                        i_event_type = 2;
                    end
                end

            end

            display(['Event Type:  ' event_type{i_event_type}]);

            %------------ try to find event type
            %h = figure(1);
            %semilogy(ot.t/year,ot.v);
            %hold on
            %semilogy(ot.t/year,ones(size(ot.t/year))*Vdyn,'r--');
            %semilogy(ot.t/year,ones(size(ot.t/year))*p.V_SS,'--','color',[0.6 0.6 0.6]);
            %xlabel('Time: (Years)')
            %ylabel('Vmax: (m/s)')
            %title([filename '   Type: ' event_type{i_event_type}],'Interpreter','none');
            %print(h,'-depsc2',filename_eps);
            %xlim([twm*(1.0 - t_end_p) twm]);
            %print(h,'-depsc2',filename_eps_z);
            %clf
            %
            L_rup_max = 0;
            R_rup_max = 0;
            Len_rup_max = 0;
            Dc_rup_max = 0;
            Lc_rup_max = 0;
            L_rup_min = 0;
            R_rup_min = 0;
            Len_rup_min = 0;
            Dc_rup_min = 0;
            Lc_rup_min = 0;           
            L_rup_mean = 0;
            R_rup_mean = 0;
            Len_rup_mean = 0;
            Dc_rup_mean = 0;
            Lc_rup_mean = 0;    
            
            if i_event_type == 5 || i_event_type ==6
                disp('Recording rupture lengths and Dc in rupture area');
                sL_rup = zeros(size(pks));
                sR_rup = zeros(size(pks));
                sLen_rup = zeros(size(pks));
                sDc_rup = zeros(size(pks));
                sLc_rup = zeros(size(pks)); 
                sVmax = zeros(size(pks));
                sT = zeros(size(pks));

                for iipks =  1:1:numel(pks)
                    sVmax(iipks) = pks(iipks);
                    sT(iipks) = locs(iipks);
                    if sVmax(iipks) >= v_th*Vdyn 
                        oxvmax = max(ox.v);
                        id0 = find(oxvmax(ox.t<= sT(iipks)) <= v_th*Vdyn,1,'last');
                        id1 = find(oxvmax(ox.t>= sT(iipks)) <= v_th*Vdyn,1,'first');
                        id1 = id1 + find((ox.t>= sT(iipks)),1,'first');
                        if id1 > numel(oxvmax)
                            id1 = numel(oxvmax);
                        end
                        ttvmax = max(ox.v(:,id0:id1),[],2);
                        iXL = find(ttvmax >= v_th*Vdyn,1,'first');
                        iXR = find(ttvmax >= v_th*Vdyn,1,'last');
                        sL_rup(iipks) = p2.X(iXL);
                        sR_rup(iipks) = p2.X(iXR);
                        sLen_rup(iipks) = sR_rup(iipks) - sL_rup(iipks); 
                        sDc_rup(iipks) = mean(p2.DC(iXL:iXR));
                        sLc_rup(iipks) = p2.MU*sDc_rup(iipks)/(p2.SIGMA*(p2.B-p2.A));                        
                    end
                   
                % find largest event
                [max_len,imax_len] = max(sLen_rup);
                L_rup_max = sL_rup(imax_len);
                R_rup_max = sR_rup(imax_len);
                Len_rup_max = sLen_rup(imax_len);
                Dc_rup_max = sDc_rup(imax_len);
                Lc_rup_max = sLc_rup(imax_len);
               
                [min_len,imin_len] = min(sLen_rup);
                L_rup_min = sL_rup(imin_len);
                R_rup_min = sR_rup(imin_len);
                Len_rup_min = sLen_rup(imin_len);
                Dc_rup_min = sDc_rup(imin_len);
                Lc_rup_min = sLc_rup(imin_len); 

                L_rup_mean = mean(sL_rup);
                R_rup_mean = mean(sR_rup);
                Len_rup_mean = mean(sLen_rup);
                Dc_rup_mean = mean(sDc_rup);
                Lc_rup_mean = mean(sLc_rup); 

                if L_rup_min < L_rup_max * 0.8
                   i_event_type = 5;
                end                                
                
		end
            end
            
            
            save(filename)

            system(['scp ' filename ' ' tardir]);
            system(['rm ' filename]);
            display([filename ' has been copyed to ' tardir]);
            %system(['scp ' filename_eps ' ' tardir]);
            %system(['rm ' filename_eps]);
            %display([filename_eps ' has been copyed to ' tardir]);
            %system(['scp ' filename_eps_z ' ' tardir]);
            %system(['rm ' filename_eps_z]);
            %display([filename_eps_z ' has been copyed to ' tardir]);  

            fid=fopen(filename_sum,'a');
            fprintf(fid,'%10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10d\n',...
                dcsigma,DC_mean,col_l,L_rup_max,R_rup_max,Len_rup_max,Dc_rup_max,Lc_rup_max,...
                L_rup_min,R_rup_min,Len_rup_min,Dc_rup_min,Lc_rup_min,...
                L_rup_mean,R_rup_mean,Len_rup_mean,Dc_rup_mean,Lc_rup_mean,vend_max,i_event_type);
            fclose(fid);

            dia.dcsigma(ii_count_glb) = dcsigma;
            dia.DC_mean(ii_count_glb) = DC_mean;
            dia.col_l(ii_count_glb) = col_l;

            dia.L_rup_max(ii_count_glb) = L_rup_max;
            dia.R_rup_max(ii_count_glb) = R_rup_max;
            dia.Len_rup_max(ii_count_glb) = Len_rup_max;
            dia.Dc_rup_max(ii_count_glb) =  Dc_rup_max;
            dia.Lc_rup_max(ii_count_glb) =  Lc_rup_max;
            dia.L_rup_min(ii_count_glb) = L_rup_min;
            dia.R_rup_min(ii_count_glb) = R_rup_min;
            dia.Len_rup_min(ii_count_glb) = Len_rup_min;
            dia.Dc_rup_min(ii_count_glb) =  Dc_rup_min;
            dia.Lc_rup_min(ii_count_glb) =  Lc_rup_min;                                   
            dia.L_rup_mean(ii_count_glb) = L_rup_mean;
            dia.R_rup_mean(ii_count_glb) = R_rup_mean;
            dia.Len_rup_mean(ii_count_glb) = Len_rup_mean;
            dia.Dc_rup_mean(ii_count_glb) =  Dc_rup_mean;
            dia.Lc_rup_mean(ii_count_glb) =  Lc_rup_mean;                                   
                        
            dia.vmax(ii_count_glb) = vend_max;
            dia.event_type(ii_count_glb) = i_event_type;
            
        end
            
            
    end

            

            clear ot ox
              
end



 

        system(['scp ' filename_sum ' ' tardir]);
        save([filename_sum '.mat']);
        system(['scp ' filename_sum '.mat ' tardir]);




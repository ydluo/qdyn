% Test of Dc slope controlling rupture

clear;
clc;

%------------------------------
rand(1,floor(sum(100*clock)));
%------------------------------


L = 100e3;
N = 2048;
DCC = 0.01;
B = 0.01;
A = 0.008;
FINITE = 1;

twm = 4000;     %simu time in years 

DCmax = 100;    % max DC cap
%HC = [4:1:10 0.2:1:5 0.4:1:5 0.6:1:5 0.8:1:5];   % L/Lc in the center
%L_slp = [0.1:0.1:1.0 1.2:0.2:2.0 2.4:0.4:4.0 5:1:10 15:5:50 60:20:200 300:100:1000];        % slope left dDc/dx normalized by critical DDCC = sigma(b-a)/mu

HC = [3];   % L/Lc in the center
L_slp = [1.0 1.2:0.2:2.0 2.4:0.4:4.0 5:1:10 15:5:50 60:20:200 300:100:1000];        % slope left dDc/dx normalized by critical DDCC = sigma(b-a)/mu
% HC = [1];   % L/Lc in the center
% L_slp = [0.3];        % slope left dDc/dx normalized by critical DDCC = sigma(b-a)/mu


%-----------


tardir = ' luoyd@pacha:/export/raid1/luoyd/qdyn_results_DC_slope/New2/';
filename_sum = 'Dc_slope_finiti_New2.txt';
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


for iiHC = 1:1:numel(HC)
    tHC = HC(iiHC);
    
    for iiL_slp = 1:1:numel(L_slp)
        tL_slp = L_slp(iiL_slp);
        
        tR_slp = tL_slp;  
            
           
            ii_count_glb = ii_count_glb + 1;
            
            filename = ['Dc_slope_finite_Lslp' num2str(tL_slp) '_Rslp' num2str(tR_slp) ...
                '_cLc' num2str(tHC) '.mat'];
            
            disp(filename);

            year = 3600*24*365;
            %if exist('qdyn')~=2, addpath ~/2D_RUPTURE/RATE_AND_STATE/qdyn/ ; end
            p = qdyn('set');
            p.B = B;
            p.A = A;
            AB_RATIO = p.A/p.B;
            p.DC = DCC;
            p.W = 100e3;

            Lb = p.MU*DCC/p.SIGMA/p.B;
            Lnuc = 1.3774*Lb;
            Lc = Lb/(1-AB_RATIO);
            Linf = 1/pi *(1-AB_RATIO)^2 *Lb;

            p.L = L;
            p.FINITE=FINITE;
            p.N = N;

            dx=p.L/p.N;
            Lb_over_dx = Lb/dx
            p.ACC = 1e-10;

            p = qdyn('set',p);


            p.TMAX = twm * year; % 
            p.NTOUT=100;
            p.NXOUT=1;
            p.NSTOP=0;
            
            p.DC = ones(size(p.X))*p.DC;
            Hcc = tHC * Lc;
            DDCC = p.SIGMA *(p.B-p.A)/p.MU * 2;
                        
            p.DC(p.X >= Hcc/2) = p.DC(p.X >= Hcc/2) + DDCC*(p.X(p.X >= Hcc/2) - Hcc/2)*tR_slp;
            p.DC(p.X <= -Hcc/2) = p.DC(p.X <= -Hcc/2) - DDCC*(p.X(p.X <= -Hcc/2) + Hcc/2)*tL_slp;
            p.DC = min(p.DC,DCmax);
            
            p.V_0 = ( 1.+0.01*exp( -(p.X/Hcc).^6 ) )*p.V_SS ;
            p.V_0 = p.V_0/mean(p.V_0)*p.V_SS;
            [p,ot,ox] = qdyn('run',p) ;
            
            Vdyn=2*mean(p.A.*p.SIGMA./p.MU.*p.VS);

            filename_eps = [filename '.eps'];
            filename_eps_z = [filename '_zoom.eps'];

            
            
                  %------------ try to find event type

            ii_end = find(ot.t >= twm*year*(1.0 - t_end_p) , 1);
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

            if vend_max <= p.V_SS*v_th_ss
                    i_event_type = 2;   %steady slip
            end

            if (vend_max > p.V_SS*v_th_ss) && (vend_max < v_th * Vdyn) || (isp == 1)
                    i_event_type = 4;   % chara aseis
                [pks,locs] = findpeaks(vend,tend,'MinPeakDistance',t_pk_dist,'MinPeakHeight',vend_max * 0.01);
                t_rec = diff(locs);
                if min(t_rec) < max(t_rec) * (1.0 - t_rec_th)   %irregular aseismic
                    i_event_type = 3;
                end

                if (vend_max <= p.V_SS*100)
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

                    if (vend_max <= p.V_SS*2)
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
                sL_rup = 0;
                sR_rup = 0;
                sLen_rup = 0;
                sDc_rup = 0;
                sLc_rup = 0; 
                iipks_2 = 0;
                sVmax = zeros(size(pks));
                sT = zeros(size(pks));

                for iipks =  1:1:numel(pks)
                    sVmax(iipks) = pks(iipks);
                    sT(iipks) = locs(iipks);
                    if sVmax(iipks) >= v_th*Vdyn 
                        iipks_2 = 1+ iipks_2;
                        oxvmax = max(ox.v);
                        id0 = find(oxvmax(ox.t<= sT(iipks)) <= v_th*Vdyn,1,'last');
                        id1 = find(oxvmax(ox.t>= sT(iipks)) <= v_th*Vdyn,1,'first');
                        id1 = id1 + find((ox.t>= sT(iipks)),1,'first');
                        if id1 > numel(oxvmax)
                            id1 = numel(oxvmax);
                        end
                        
                        ttvmax = max(ox.v(:,id0:id1),[],2);
                        if max(ttvmax) >= v_th*Vdyn 
                            iipks_2 = 1+ iipks_2;
                            iXL = find(ttvmax >= v_th*Vdyn,1,'first');
                            iXR = find(ttvmax >= v_th*Vdyn,1,'last');
                            sL_rup(iipks_2) = p.X(iXL);
                            sR_rup(iipks_2) = p.X(iXR);
                            sLen_rup(iipks_2) = sR_rup(iipks_2) - sL_rup(iipks_2); 
                            sDc_rup(iipks_2) = mean(p.DC(iXL:iXR));
                            sLc_rup(iipks_2) = p.MU*sDc_rup(iipks_2)/(p.SIGMA*(p.B-p.A));  
                        end
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
               	
		if Len_rup_min < Len_rup_max * 0.8
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
                tHC,tL_slp,tR_slp,L_rup_max,R_rup_max,Len_rup_max,Dc_rup_max,Lc_rup_max,...
                L_rup_min,R_rup_min,Len_rup_min,Dc_rup_min,Lc_rup_min,...
                L_rup_mean,R_rup_mean,Len_rup_mean,Dc_rup_mean,Lc_rup_mean,vend_max,i_event_type);
            fclose(fid);

            dia.HC(ii_count_glb) = tHC;
            dia.L_slp(ii_count_glb) = tL_slp;
            dia.R_slp(ii_count_glb) = tR_slp;
            
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

            
            disp(['Done'])
          %  clear ot ox
              
end



 

        system(['scp ' filename_sum ' ' tardir]);
        save([filename_sum '.mat']);
        system(['scp ' filename_sum '.mat ' tardir]);




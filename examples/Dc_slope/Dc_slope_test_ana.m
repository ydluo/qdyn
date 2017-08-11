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
HC_new = [0:1:10 0.2:1:5 0.4:1:5 0.6:1:5 0.8:1:5];   % L/Lc in the center
L_slp_new = [0.1:0.1:1.0 1.2:0.2:2.0 2.4:0.4:4.0 5:1:10 15:5:50 60:20:200 300:100:1000];        % slope left dDc/dx normalized by critical DDCC = sigma(b-a)/mu

% HC = [1];   % L/Lc in the center
% L_slp = [0.3];        % slope left dDc/dx normalized by critical DDCC = sigma(b-a)/mu


tardir_new = ' luoyd@pacha:/export/raid1/luoyd/qdyn_results_DC_slope/New2/';
filename_sum_new = 'New2_Dc_slope_finite.txt';
event_type_new = {'Undefined','Steady-Slip','Irreg. Aseismic','Chara. Aseismic','Irreg. Seismic','Chara. Seismic'};

t_pk_dist_new = 1800;        % min peak distance in seconds

v_th_new = 0.1; %threshold for dynamic events;
v_th_ss_new = 1.01; 
t_rec_th_new = 0.2; %threshold for irregular events; t_rec_min < t_rec_max * (1-t_rec_th)
num_cc_new = 5;     % number of bins used to decide if a event is cont. dec.
t_end_p_new = 0.5;  %portion of twm used for analysiscd 


ii_count_glb_new = 0;
% dia.HC = zeros(ttl_sim,1);
% dia.L_slp = zeros(ttl_sim,1);
% dia.R_slp = zeros(ttl_sim,1);
% dia.L_rup = zeros(ttl_sim,1);
% dia.R_rup = zeros(ttl_sim,1);
% dia.Dc_rup = zeros(ttl_sim,1);
% dia.Lc_rup = zeros(ttl_sim,1);
% dia.vmax = zeros(ttl_sim,1);
% dia.event_type = zeros(ttl_sim,1);


for iiHC = 1:1:numel(HC_new)
    tHC = HC_new(iiHC);
    
    for iiL_slp = 1:1:numel(L_slp_new)
        tL_slp = L_slp_new(iiL_slp);
        
        tR_slp = tL_slp;  
            
           
            ii_count_glb_new = ii_count_glb_new + 1;
            
            ffname = ['Dc_slope_finite_Lslp' num2str(tL_slp) '_Rslp' num2str(tR_slp) ...
                '_cLc' num2str(tHC) '.mat'];
            
            disp(['Processing' ffname ' ...']);

            load(ffname);
            
            Vdyn=2*mean(p.A.*p.SIGMA./p.MU.*p.VS);

            filename_eps = [ffname '.eps'];
            filename_eps_z = [ffname '_zoom.eps'];
            filename_rup = [ffname '_rup.eps'];

            
             %------------ try to find event type

            
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
                i_event_type = 6;
                sL_rup = 0;
                sR_rup = 0;
                sLen_rup = 0;
                sDc_rup = 0;
                sLc_rup = 0; 
                iipks_2 = 0;
                sVmax = zeros(size(pks));
                sT = zeros(size(pks));
                sT2 = 0;

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
                        if max(ttvmax) >= v_th*Vdyn 
                            iipks_2 = 1+ iipks_2;
                            sT2(iipks_2) = sT(iipks);
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

                if Len_rup_min < Len_rup_max * 0.95
                   i_event_type = 5;
                end
                
                end

                h2 = figure(2);
                for ii_pt = 1:1:numel(sLen_rup)
                    plot([sT2(ii_pt) sT2(ii_pt)]/year,[sL_rup(ii_pt) sR_rup(ii_pt)]/1000,'r','Linewidth',1);
                    hold on
                    plot([sT2(ii_pt) sT2(ii_pt)]/year,([-.5*sLc_rup(ii_pt) .5*sLc_rup(ii_pt)]+(sL_rup(ii_pt)+sR_rup(ii_pt))/2)/1000,...
                        'b','Linewidth',1);  
                end
    %             semilogy(ot.t/year,ones(size(ot.t/year))*Vdyn,'r--');
    %             semilogy(ot.t/year,ones(size(ot.t/year))*p.V_SS,'--','color',[0.6 0.6 0.6]);
                xlabel('Time: (Years)')
                ylabel('X: (km)')
                title([ffname '   Type: ' event_type{i_event_type}],'Interpreter','none');
                xlim([twm*(1.0 - t_end_p) twm]);
                ylim([min(p.X) max(p.X)]/1000);
                legend('Rupture Area','Lc');
                print(h2,'-depsc2',filename_rup);

                clf                  
                
            end
            
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
            
            h = figure(1);
            semilogy(ot.t/year,ot.v);
            hold on
            semilogy(ot.t/year,ones(size(ot.t/year))*Vdyn,'r--');
            semilogy(ot.t/year,ones(size(ot.t/year))*p.V_SS,'--','color',[0.6 0.6 0.6]);
            xlabel('Time: (Years)')
            ylabel('Vmax: (m/s)')
            title([ffname '   Type: ' event_type{i_event_type}],'Interpreter','none');
            print(h,'-depsc2',filename_eps);
            xlim([twm*(1.0 - t_end_p) twm]);
            print(h,'-depsc2',filename_eps_z);
            clf

   
            
            clear ot ox
            
    end
              

    

end
 

        system(['scp ' filename_sum_new ' ' tardir_new]);
        save([filename_sum_new '.mat']);
        system(['scp ' filename_sum_new '.mat ' tardir_new]);




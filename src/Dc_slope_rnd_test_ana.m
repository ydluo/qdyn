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
bin_col_l = [250:125:1000 1250:250:5000 6e3:1e3:10e3 12e3:2e3:30e3];	%r of col 

%bin_DC_mean = [0.025];		%mean DC
%bin_dcsigma = [0.25];		%Dv
%bin_col_l = [5000];	%r of col 



DC_min = 0.01;		%min_DC




year = 3600*24*365;      
      
      
%-----------


tardir = ' luoyd@pacha:/export/raid1/luoyd/qdyn_results_DC_test_rnd/';
filename_sum_new = 'New_Dc_test_rnd_finiti.txt';
event_type = {'Undefined','Steady-Slip','Irreg. Aseismic','Chara. Aseismic','Irreg. Seismic','Chara. Seismic'};

t_pk_dist = 1800;        % min peak distance in seconds

v_th = 0.1; %threshold for dynamic events;
v_th_ss = 1.01; 
t_rec_th = 0.2; %threshold for irregular events; t_rec_min < t_rec_max * (1-t_rec_th)
num_cc = 5;     % number of bins used to decide if a event is cont. dec.
t_end_p = 0.5;  %portion of twm used for analysiscd 

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



            
for iibin_dcsigma = 1:1:numel(bin_dcsigma)
    dcsigma = bin_dcsigma(iibin_dcsigma);
    
    for iiDC_mean = 1:1:numel(bin_DC_mean)
        DC_mean = bin_DC_mean(iiDC_mean);
        
        for iicol_l = 1:1:numel(bin_col_l)
            col_l = bin_col_l(iicol_l);
                    
           
            ii_count_glb_new = ii_count_glb_new + 1;
            
            DC_v = (exp((dcsigma)^2)-1)*DC_mean^2;
            dcmu = log(DC_mean^2/sqrt(DC_v+DC_mean^2));
            
            filename = ['Dc_test_rnd_finite_dcsigma' num2str(dcsigma) '_DCmean' num2str(DC_mean) ...
                '_coll' num2str(col_l) '.mat'];
            
            disp(['Processing' filename ' ...']);

            load(filename);
            
            Vdyn=2*mean(p.A.*p.SIGMA./p.MU.*p.VS);

            filename_eps = [filename '.eps'];
            filename_eps_z = [filename '_zoom.eps'];           
            

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
                            sL_rup(iipks_2) = p2.X(iXL);
                            sR_rup(iipks_2) = p2.X(iXR);
                            sLen_rup(iipks_2) = sR_rup(iipks_2) - sL_rup(iipks_2); 
                            sDc_rup(iipks_2) = mean(p2.DC(iXL:iXR));
                            sLc_rup(iipks_2) = p2.MU*sDc_rup(iipks_2)/(p2.SIGMA*(p2.B-p2.A)); 
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
            
            

            fid=fopen(filename_sum_new,'a');
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
 
            h = figure(1);
            semilogy(ot.t/year,ot.v);
            hold on
            semilogy(ot.t/year,ones(size(ot.t/year))*Vdyn,'r--');
            semilogy(ot.t/year,ones(size(ot.t/year))*p.V_SS,'--','color',[0.6 0.6 0.6]);
            xlabel('Time: (Years)')
            ylabel('Vmax: (m/s)')
            title([filename '   Type: ' event_type{i_event_type}],'Interpreter','none');
            print(h,'-depsc2',filename_eps);
            xlim([twm2*(1.0 - t_end_p) twm2]);
            print(h,'-depsc2',filename_eps_z);
            clf
                               
            
        end
            
            
    end


   

            clear ot ox
              
end



 

        system(['scp ' filename_sum_new ' ' tardir]);
        save([filename_sum_new '.mat']);
        system(['scp ' filename_sum_new '.mat ' tardir]);
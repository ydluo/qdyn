label = {'000', '025', '050', '075', '100', '125', '150', '175', '200', '250', '300', '350'};
NT = 1;
for h = 1:length(label)
    fid = fopen([ 'fltst_dp', label{h} ], 'w');
    fprintf(fid, '# problem=SEAS Benchmark No.1\n');
    fprintf(fid, '# author=Y.Luo, B.Idini, J.-P.Ampuero\n');
    fprintf(fid, '# date=2018/04/19\n');
    fprintf(fid, '# code=qdyn\n');
    fprintf(fid, '# code_version=1.1\n');
    fprintf(fid, ['# element_size=', num2str(p.L/p.N), ' (m)\n'] );
    fprintf(fid, ['# minimum_time_step=', num2str(min(diff(ot.t(1:NT:end)))), ' s\n'] );
    fprintf(fid, ['# maximum_time_step=', num2str(max(diff(ot.t(1:NT:end)))/year), ' yr\n'] );
    fprintf(fid, ['# num_time_steps=', num2str(length(ot.t)), '\n'] );
    fprintf(fid, ['# location= on fault, ', num2str(depth(h)), 'km depth\n']);
    fprintf(fid, '# Column #1 = Time (s)\n');
    fprintf(fid, '# Column #2 = Slip (m)\n');
    fprintf(fid, '# Column #3 = Slip rate (log10 m/s)\n');
    fprintf(fid, '# Column #4 = Shear stress (MPa)\n');
    fprintf(fid, '# Column #5 = State (log10 s)\n');
    fprintf(fid, '#\n');
    fprintf(fid, 't slip slip_rate shear_stress state\n');
 
    fid0 = fopen(['fort.10', num2str(h,'%03d')], 'r');
    fgetl(fid0);
    fgetl(fid0);
    rdat = textscan(fid0, '%f %f %f %f %f %f');
    fclose(fid0);
    t = rdat{1};
    V = log10(rdat{2});
    s = log10(rdat{3});
    tau = rdat{4}/1e6;
    D = rdat{5};
    fprintf(fid, '%.15e %.15e %.15e %.15e %.15e %.15e\n', [t(1:NT:end); D(1:NT:end); V(1:NT:end); tau(1:NT:end); s(1:NT:end)]); 
    figure; plot(t(1:NT:end), tau(1:NT:end))
    fclose(fid);
end

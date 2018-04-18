label = {'000', '025', '050', '075', '100', '125', '150', '175', '200', '250', '300', '350'};
for h = 1:9
    fid = fopen([ 'fltst_dp', label{h} ], 'w');
    fprintf(fid, '# problem=SEAS Benchmark No.1\n');
    fprintf(fid, '# author=Y.Luo, B.Idini, J.-P.Ampuero\n');
    fprintf(fid, '# date=2018/04/19\n');
    fprintf(fid, '# code=qdyn\n');
    fprintf(fid, '# code_version=1.1\n');
    fprintf(fid, ['# element_size=\n', num2str(p.L/p.N), ' (m)'] );
    fprintf(fid, ['# minimum_time_step=\n', num2str(min(diff(ot.t))] );
    fprintf(fid, ['# maximum_time_step=\n', num2str(max(diff(ot.t))] );
    fprintf(fid, ['# num_time_steps=\n', num2str(length(ot.t))] );
    fprintf(fid, ['# location= on fault, ', num2str(depth(t)), 'km depth\n']);
    fprintf(fid, 'Column #1 = Time (s)\n');
    fprintf(fid, 'Column #2 = Slip (m)\n');
    fprintf(fid, 'Column #3 = Slip rate (log10 m/s)\n');
    fprintf(fid, 'Column #4 = Shear stress (MPa)\n');
    fprintf(fid, 'Column #5 = State (log10 s)\n');
    fprintf(fid, '#\n');
    fprintf(fid, 't slip slip_rate shear_stress state\n');
 
    fid0 = open(['fort.1000', num2str(h)], 'r');
    rdat = textscan(fid0, '%f %f %f %f %f %f'); 
    t = rdat{1};
    V = log10(rdat{2});
    s = log10(rdat{3});
    tau = rdat{4}/1e6;
    D = rdat{5};
    fprintf(fid, '%.15g %.15g %.15g %.15g %.15g %.15g\n', [t; D; V; tau; s]); 
    fclose(fid);
end

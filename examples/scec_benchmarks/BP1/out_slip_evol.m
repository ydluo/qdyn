fid = fopen('devol.dat', 'w'); 
fprintf(fid, '# problem=SEAS Benchmark No.1\n');
fprintf(fid, '# author=Y.Luo, B.Idini, J.-P.Ampuero\n');
fprintf(fid, '# date=2018/04/19\n');
fprintf(fid, '# code=qdyn\n');
fprintf(fid, '# code_version=1.1\n');
fprintf(fid, ['# element_size=', num2str(p.L/p.N),' (m)\n'] );
fprintf(fid, '# Row #1 = Depth (m) with two zeros first\n');
fprintf(fid, '# Column #1 = Time (s)\n');
fprintf(fid, '# Column #2 = Max slip rate (m/s)\n');
fprintf(fid, '# Column #3-83 = Slip (m)\n');
fprintf(fid, '#\n');
fprintf(fid, 'z\n');
fprintf(fid, 't max_slip_rate slip\n');
 
t = ox.t;
z = ox.x + 20000;
v_max = max(ox.v, [], 1);
D = ox.d;

fprintf(fid, '%.15g\t', [0, 0, z']);
fprintf(fid, '\n');
for i = 1:length(t)
    fprintf(fid, '%.15e\t', [t(i), v_max(i), D(:, i)']);
    fprintf(fid, '\n');
end
fclose(fid);


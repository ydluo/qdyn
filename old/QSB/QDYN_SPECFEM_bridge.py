import os

n_core_allco = 108;
n_core_allco_single_node = 12;

QDYN_dir_work = "/home/luoyd/qdyn_git_svn/trunk/src/"
QDYN_dir_out_store = "/home/luoyd/QSB_out/test/QDYN_out/"
SPECFEM_dir_work = "/home/luoyd/SPECFEM3D_CIG_RSF_10/bin/"
SPECFEM_dir_in = "/home/luoyd/SPECFEM3D_CIG_RSF_10/DATA/"
SPECFEM_dir_out = "/home/luoyd/SPECFEM3D_CIG_RSF_10/OUTPUT_FILES/"
SPECFEM_dir_out_store = "/home/luoyd/QSB_out/test/SPECFEM_out/"

cmd = "cp SEM_to_QDYN_RSF.m "+SPECFEM_dir_out 
os.system(cmd)
cmd = "cp FSEM3D_snapshot.m "+SPECFEM_dir_out  
os.system(cmd)
cmd = "cp Qdyn_read_ox_seq.m "+SPECFEM_dir_out
os.system(cmd)
cmd = "cp Qdyn_read_in.m "+SPECFEM_dir_out
os.system(cmd)
cmd = "cp nodesonfault "+QDYN_dir_work
os.system(cmd)

for irun in range(1,5):

 os.system("date")
 print "QSB: run # %d ..." % (irun)
 cmd = "cd "+ QDYN_dir_work
 os.system(cmd)
 print "QSB: run # %d QDYN simulation ..." % (irun)
 cmd = "mpirun -np "+str(n_core_allco_single_node)+ " qdyn"
 os.system(cmd)
 print "QSB: QDYN run # %d finished" % (irun)
 cmd = "cp qdyn.in "+SPECFEM_dir_out
 os.system(cmd)
 cmd = "cp fort.20001 "+SPECFEM_dir_out
 os.system(cmd)

 os.system("date")
 print "QSB: Matlab run # %d QDYN to SEM ..." % (irun)
 os.system("matlab -nosplash -nodesktop -r 'QDYN_to_SEM_m_RSF_f; quit'")
 print "QSB: Matlab run # %d QDYN to SEM finished" % (irun)
 cmd = "cp ./input_file.txt "+SPECFEM_dir_in
 os.system(cmd)
 print "QSB: input_file.txt moved to "+SPECFEM_dir_in

 os.system("date")
 ndir = QDYN_dir_out_store+"run"+str(irun) 
 cmd = "mkdir "+ndir
 os.system(cmd)
 cmd = "mv " +QDYN_dir_work+ " fort.* "+ndir
 os.system(cmd)
 cmd = "cp " +QDYN_dir_work+ " timpestamp.txt "+ndir
 os.system(cmd)
 print "QSB: QDYN outputs moved to "+ndir

 os.system("date")
 print "QSB: run # %d SPECFEM simulation ..." % (irun)
 cmd = "cd "+SPECFEM_dir_work
 os.system(cmd)
 cmd = "mpirun -np "+str(n_core_allco)+" ./xspecfem3D"
 os.system(cmd)
 os.system("date")
 print "QSB: SPECFEM run # %d finished" % (irun)
 
 cmd = "cd "+SPECFEM_dir_out
 os.system("date")
 print "QSB: Matlab run # %d SEM to QDYN ..." % (irun)
 os.system(cmd)
 os.system("matlab -nosplash -nodesktop -r 'SEM_to_QDYN_RSF; quit'")
 print "QSB: Matlab run # %d SEM to QDYN finished" % (irun)
 cmd = "cp qdyn.in "+QDYN_dir_work
 os.system(cmd)

 os.system("date")
 ndir = SPECFEM_dir_out_store+"run"+str(irun)
 cmd = "mkdir "+ndir
 os.system(cmd)
 cmd = "mv " +SPECFEM_dir_out_store+ "/* "+ndir
 os.system(cmd)
 print "QSB: SPECFEM outputs moved to "+ndir



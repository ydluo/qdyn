#!/bin/bash

N_core_allco=108
N_core_allco_single_node=12

QDYN_dir_work="/home/luoyd/qdyn_git_svn/trunk/src/"
QDYN_dir_out_store="/home/luoyd/QSB_out/test/QDYN_out/"
SPECFEM_dir_work="/home/luoyd/SPECFEM3D_CIG_RSF_10/bin/"
SPECFEM_dir_in="/home/luoyd/SPECFEM3D_CIG_RSF_10/DATA/"
SPECFEM_dir_out="/home/luoyd/SPECFEM3D_CIG_RSF_10/OUTPUT_FILES/"
SPECFEM_dir_out_store="/home/luoyd/QSB_out/test/SPECFEM_out/"

cp SEM_to_QDYN_RSF.m $SPECFEM_dir_ou &&
cp FSEM3D_snapshot.m $SPECFEM_dir_out &&
cp Qdyn_read_ox_seq.m $SPECFEM_dir_out &&
cp cp Qdyn_read_in.m $SPECFEM_dir_out &&
cp nodesonfault $QDYN_dir_work &&


echo "N_core_allco_single_node = $N_core_allco_single_node &&
echo "N_core_allco = $N_core_allco" &&


for irun in {1..2}
do

   date &&
   echo QSB: run no. $irun &&
   cd $QDYN_dir_work &&
   echo QSB: run no. $irun QDYN simulation ... &&
   ./qdyn &&
   echo QSB: run no. $irun QDYN simulation finished &&
   cp qdyn.in $SPECFEM_dir_out &&
   cp fort.20001 $SPECFEM_dir_out &&
   
   date &&
   echo QSB: Matlab run no. $irun QDYN to SEM ... &&
   matlab -nosplash -nodesktop -r 'QDYN_to_SEM_m_RSF_f; quit' &&
   echo QSB: Matlab run no. $irun QDYN to SEM finished &&
   cp ./input_file.txt $SPECFEM_dir_in &&
   echo QSB: input_file.txt moved to $SPECFEM_dir_in &&
   
   date &&
   mkdir ${QDYN_dir_out_store}/run$irun &&
   mv ${QDYN_dir_work}/fort.* ${QDYN_dir_out_store}/run$irun &&
   cp ${QDYN_dir_work}/timpestamp.txt ${QDYN_dir_out_store}/run$irun &&
   echo QSB: QDYN outputs moved to ${QDYN_dir_out_store}/run$irun &&
   
   date &&
   echo "QSB: run no. $irun SPECFEM simulation ... &&
   cd $SPECFEM_dir_work &&
   mpirun -np $N_core_allco ./xspecfem3D &&
   echo "QSB: run no. $irun SPECFEM simulation finished &&
   
   cd $SPECFEM_dir_out &&
   date &&
   echo QSB: Matlab run no. $irun SEM to QDYN ... &&
   matlab -nosplash -nodesktop -r 'SEM_to_QDYN_RSF; quit' &&
   cp qdyn.in $QDYN_dir_work &&
   echo QSB: Matlab run no. $irun SEM to QDYN finished &&

   date &&
   mkdir ${SPECFEM_dir_out_store}/run$irun &&
   mv ${SPECFEM_dir_out}/* ${SPECFEM_dir_out_store}/run$irun &&
   echo QSB: SPECFEM outputs moved to ${SPECFEM_dir_out_store}/run$irun

done

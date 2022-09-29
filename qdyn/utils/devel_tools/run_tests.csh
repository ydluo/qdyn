#!/bin/csh -f
# runs a test in each directory containing a matlab script named analyze_test.m

echo
echo "-----------------------------"
echo "---  Test suite for QDYN  ---"
echo "-----------------------------"
echo

@ ntotal = 0
@ nok = 0
set start=`pwd`
foreach name (`ls -d test_*`)
  if (! -d $name) continue
  if (! -e $name/run_and_analyze_test.m) continue

  @ ntotal ++
  cd $name 
  echo -n $name "............"
    
 # run test 
  rm -f test.out
  matlab -nodisplay -nojvm -nosplash -nodesktop < run_and_analyze_test.m > test.out

  grep -sq "Test = 1" test.out  
 # $status is exit code from last command 
 # $status==0 means succesful
  if ($status == 0) then
    @ nok ++
    echo "......... [OK]"
  else
    echo "......... [FAILED]"
  endif

  cd $start
end

echo
echo "Summary: "
if ($nok == $ntotal) then
  echo "[OK] Passed all " $ntotal " tests"
else
  echo "[FAILED] Only passed " $nok " tests out of " $ntotal 
endif

echo

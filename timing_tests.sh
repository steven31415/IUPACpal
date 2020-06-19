#!/bin/bash

start_measuring_time() {
  read s1 s2 < <(date +'%s %N')
}

stop_measuring_time() {
  read e1 e2 < <(date +'%s %N')
}

get_elapsed_time() {
  TIME_TAKEN=$(bc <<< "scale=3; $((e1-s1))+($((10#$e2-10#$s2))/1000000000)")
  result=$TIME_TAKEN
}

today=`date +%Y-%m-%d.%H:%M:%S`
printf "test_no\tfile\tseq_name\tmin_length\tmax_length\tmax_gap\tmismatches\tIUPAC_correct_pals\tIUPAC_incorrect_pals\tIUPAC_runtime\n" > "test_results_$today.csv"

mkdir temp
FILENAME="timing_tests.cfg"
TEST_COUNT=1

while read -r LINE
do
  if [[ $LINE == "file "* ]]; then
  input_file=${LINE:5}
  fi

  if [[ $LINE == "s "* ]]; then
  s=${LINE:2}
  fi

  if [[ $LINE == "m "* ]]; then
  m=${LINE:2}
  fi

  if [[ $LINE == "M "* ]]; then
  M=${LINE:2}
  fi

  if [[ $LINE == "g "* ]]; then
  g=${LINE:2}
  fi

  if [[ $LINE == "x "* ]]; then
  x=${LINE:2}
  fi

  if [[ $LINE == "" ]]; then
    echo "RUNNING TEST NO: $TEST_COUNT"

    start_measuring_time
    ./IUPACpal -f $input_file -s $s -m $m -M $M -g $g -x $x -o temp/IUPACpal.out &> /dev/null
    stop_measuring_time
    get_elapsed_time
    IUPACPAL_TIME_TAKEN=$result
    echo "IUPACpal:              $IUPACPAL_TIME_TAKEN seconds"

    python tools/check_correctness.py $input_file $s temp/IUPACpal.out $m $M $g $x > temp/IUPAC_correctness.out
    IUPACPAL_CORRECT_PALS=$(awk '{print $2}' temp/IUPAC_correctness.out | sed -n '9p')
    IUPACPAL_INCORRECT_PALS=$(awk '{print $2}' temp/IUPAC_correctness.out | sed -n '10p')
    echo "Correct in IUPACpal:   $IUPACPAL_CORRECT_PALS"
    echo "Incorrect in IUPACpal: $IUPACPAL_INCORRECT_PALS"

    echo "Test $TEST_COUNT complete"

    printf "$TEST_COUNT\t$input_file\t$s\t$m\t$M\t$g\t$x\t$IUPACPAL_CORRECT_PALS\t$IUPACPAL_INCORRECT_PALS\t$IUPACPAL_TIME_TAKEN\n" >> "test_results_$today.csv"

    echo
    TEST_COUNT=$((TEST_COUNT+1))
  fi
done < "$FILENAME"

rm -rf temp

echo "All tests completed!"

#write output to a csv file
#include reading of palindrome output analysed by python correctness and compare tools

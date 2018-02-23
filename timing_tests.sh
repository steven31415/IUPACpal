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
printf "test_no\tfile\tseq_name\tmin_length\tmax_length\tmax_gap\tmismatches\tIUPAC_total_pals\tIUPAC_unique_pals\tIUPAC_correct_pals\tIUPAC_incorrect_pals\temboss_total_pals\temboss_unique_pals\temboss_correct_pals\temboss_incorrect_pals\tshared_pals\tIUPAC_runtime\temboss_runtime\n" > "test_results_$today.csv"

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


    start_measuring_time
    ./emboss_palindrome $input_file -minpallen $m -maxpallen $M -gaplimit $g -nummismatches $x -outfile temp/emboss_palindrome.out -overlap &> /dev/null
    stop_measuring_time
    get_elapsed_time
    EMBOSS_TIME_TAKEN=$result
    echo "emboss:                $EMBOSS_TIME_TAKEN seconds"


    python tools/compare_output.py temp/IUPACpal.out temp/emboss_palindrome.out > temp/compare.out
    TOTAL_IN_IUPACPAL=$(awk '{print $2}' temp/compare.out | sed -n '3p')
    TOTAL_IN_EMBOSS=$(awk '{print $2}' temp/compare.out | sed -n '4p')
    IN_BOTH=$(awk '{print $2}' temp/compare.out | sed -n '5p')
    ONLY_IN_IUPACPAL=$(awk '{print $2}' temp/compare.out | sed -n '6p')
    ONLY_IN_EMBOSS=$(awk '{print $2}' temp/compare.out | sed -n '7p')

    echo "Pals in IUPACpal:      $TOTAL_IN_IUPACPAL"
    echo "Pals in emboss:        $TOTAL_IN_EMBOSS"
    echo "Pals in both:          $IN_BOTH"
    echo "Only in IUPACpal:      $ONLY_IN_IUPACPAL"
    echo "Only in emboss:        $ONLY_IN_EMBOSS"


    python tools/check_correctness.py $input_file $s temp/IUPACpal.out $m $M $g $x > temp/IUPAC_correctness.out
    IUPACPAL_CORRECT_PALS=$(awk '{print $2}' temp/IUPAC_correctness.out | sed -n '9p')
    IUPACPAL_INCORRECT_PALS=$(awk '{print $2}' temp/IUPAC_correctness.out | sed -n '10p')
    echo "Correct in IUPACpal:   $IUPACPAL_CORRECT_PALS"
    echo "Incorrect in IUPACpal: $IUPACPAL_INCORRECT_PALS"

    #this correctness check assumes that f_IUPACpal and f_emboss represent that same data (albeit in different formats)
    python tools/check_correctness.py $input_file $s temp/emboss_palindrome.out $m $M $g $x > temp/emboss_correctness.out
    EMBOSS_CORRECT_PALS=$(awk '{print $2}' temp/emboss_correctness.out | sed -n '9p')
    EMBOSS_INCORRECT_PALS=$(awk '{print $2}' temp/emboss_correctness.out | sed -n '10p')
    echo "Correct in emboss:     $EMBOSS_CORRECT_PALS"
    echo "Incorrect in emboss:   $EMBOSS_INCORRECT_PALS"


    echo "Test $TEST_COUNT complete"

    printf "$TEST_COUNT\t$input_file\t$s\t$m\t$M\t$g\t$x\t$TOTAL_IN_IUPACPAL\t$ONLY_IN_IUPACPAL\t$IUPACPAL_CORRECT_PALS\t$IUPACPAL_INCORRECT_PALS\t$TOTAL_IN_EMBOSS\t$ONLY_IN_EMBOSS\t$EMBOSS_CORRECT_PALS\t$EMBOSS_INCORRECT_PALS\t$IN_BOTH\t$IUPACPAL_TIME_TAKEN\t$EMBOSS_TIME_TAKEN\n" >> "test_results_$today.csv"

    echo
    TEST_COUNT=$((TEST_COUNT+1))
  fi
done < "$FILENAME"

rm -rf temp

echo "All tests completed!"

#write output to a csv file
#include reading of palindrome output analysed by python correctness and compare tools
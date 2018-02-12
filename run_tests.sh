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

if [ $# -ne 1 ]; then 
  echo "Must provide a single filename as test configuration"
  exit -1
fi

today=`date +%Y-%m-%d.%H:%M:%S`
printf "test_no\tsoftware\tfile\tmin_length\tmax_length\tmax_gap\tmismatches\ttotal_pals\tshared_pals\tunique_pals\tcorrect_pals\tincorrect_pals\truntime\n" > "test_results_$today.csv"

FILENAME="$1"
TEST_COUNT=1
while read -r LINE
do
  if [[ $LINE == "f_emboss "* ]]; then
	f_emboss=${LINE:9}
  fi

  if [[ $LINE == "f_IUPACpal "* ]]; then
	f_IUPACpal=${LINE:11}
  fi

  if [[ $LINE == "m "* ]]; then
	m=${LINE:1}
  fi

  if [[ $LINE == "M "* ]]; then
	M=${LINE:1}
  fi

  if [[ $LINE == "g "* ]]; then
	g=${LINE:1}
  fi

  if [[ $LINE == "x "* ]]; then
	x=${LINE:1}
  fi

  if [[ $LINE == "" ]]; then
	echo "Running test $TEST_COUNT"


	start_measuring_time
	./IUPACpal -f $f_IUPACpal -m $m -M $M -g $g -x $x &> /dev/null
	stop_measuring_time
	get_elapsed_time
  IUPACPAL_TIME_TAKEN=$result
  echo "IUPACpal: $IUPACPAL_TIME_TAKEN seconds"


  start_measuring_time
  ./emboss_palindrome $f_emboss -minpallen $m -maxpallen $M -gaplimit $g -nummismatches $x -outfile emboss_palindrome.out -overlap &> /dev/null
  stop_measuring_time
  get_elapsed_time
  EMBOSS_TIME_TAKEN=$result
  echo "emboss:   $EMBOSS_TIME_TAKEN seconds"
	

  python tools/compare_output.py IUPACpal.out emboss_palindrome.out > compare.out
  TOTAL_IN_IUPACPAL=$(awk '{print $2}' compare.out | sed -n '4p')
  TOTAL_IN_EMBOSS=$(awk '{print $2}' compare.out | sed -n '5p')
  IN_BOTH=$(awk '{print $2}' compare.out | sed -n '7p')
  ONLY_IN_IUPACPAL=$(awk '{print $2}' compare.out | sed -n '8p')
  ONLY_IN_EMBOSS=$(awk '{print $2}' compare.out | sed -n '9p')

  echo "Pals in IUPACpal: $TOTAL_IN_IUPACPAL"
  echo "Pals in emboss:   $TOTAL_IN_EMBOSS"
  echo "Pals in both:     $IN_BOTH"
  echo "Only in IUPACpal: $ONLY_IN_IUPACPAL"
  echo "Only in emboss:   $ONLY_IN_EMBOSS"


  python tools/check_correctness.py $f_IUPACpal IUPACpal.out $x > correctness.out
  IUPACPAL_CORRECT_PALS=$(awk '{print $3}' correctness.out | sed -n '6p')
  IUPACPAL_INCORRECT_PALS=$(awk '{print $3}' correctness.out | sed -n '7p')
  echo "Correct in IUPACpal:   $IUPACPAL_CORRECT_PALS"
  echo "Incorrect in IUPACpal: $IUPACPAL_INCORRECT_PALS"

  #this correctness check assumes that f_IUPACpal and f_emboss represent that same data (albeit in different formats)
  python tools/check_correctness.py $f_IUPACpal emboss_palindrome.out $x > correctness.out
  EMBOSS_CORRECT_PALS=$(awk '{print $3}' correctness.out | sed -n '6p')
  EMBOSS_INCORRECT_PALS=$(awk '{print $3}' correctness.out | sed -n '7p')
  echo "Correct in emboss:   $EMBOSS_CORRECT_PALS"
  echo "Incorrect in emboss: $EMBOSS_INCORRECT_PALS"


  echo "Test $TEST_COUNT complete"

  printf "$TEST_COUNT\tIUPACpal\t$f_IUPACpal\t$m\t$M\t$g\t$x\t$TOTAL_IN_IUPACPAL\t$IN_BOTH\t$ONLY_IN_IUPACPAL\t$IUPACPAL_CORRECT_PALS\t$IUPACPAL_INCORRECT_PALS\t$IUPACPAL_TIME_TAKEN\n" >> "test_results_$today.csv"
  printf "$TEST_COUNT\temboss\t$f_emboss\t$m\t$M\t$g\t$x\t$TOTAL_IN_EMBOSS\t$IN_BOTH\t$ONLY_IN_EMBOSS\t$EMBOSS_CORRECT_PALS\t$EMBOSS_INCORRECT_PALS\t$EMBOSS_TIME_TAKEN\n" >> "test_results_$today.csv"

  echo
	echo
	TEST_COUNT=$((TEST_COUNT+1))
  fi
done < "$FILENAME"

echo "All tests completed!"

#write output to a csv file
#include reading of palindrome output analysed by python correctness and compare tools
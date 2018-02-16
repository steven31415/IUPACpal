#!/bin/bash

if [ $# -ne 1 ]; then 
  echo "Must provide a single filename to use as test input"
  exit -1
fi

FILENAME="$1"

test_values=(
        2 8 5 1
        0 0 0 0
      )

LOOP_COUNT=$(( ${#test_values[@]} / 4 - 1))

for i in `seq 0 $LOOP_COUNT`;
do
  ./IUPACpal -f $FILENAME -m ${test_values[(( 0 + i * 4 ))]} -M ${test_values[(( 1 + i * 4 ))]} -g ${test_values[(( 2 + i * 4 ))]} -x ${test_values[(( 3 + i * 4 ))]} -o test_temp1.out >/dev/null
  python tools/check_correctness.py $FILENAME test_temp1.out ${test_values[(( 0 + i * 4 ))]} ${test_values[(( 1 + i * 4 ))]} ${test_values[(( 2 + i * 4 ))]} ${test_values[(( 3 + i * 4 ))]} > test_temp2.out
  INCORRECT_COUNT=$(awk '{print $3}' test_temp2.out | sed -n '10p')
  cat test_temp2.out
  if [ $INCORRECT_COUNT = 0 ]; then
    echo "TEST $(( i + 1 )) PASSED"
  else
    echo "TEST $(( i + 1 )) FAILED"
  fi
done

rm test_temp1.out
rm test_temp2.out
echo "All tests completed!"

#write output to a csv file
#include reading of palindrome output analysed by python correctness and compare tools
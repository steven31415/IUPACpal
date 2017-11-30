#!/bin/bash
start_measuring_time() {
  read s1 s2 < <(date +'%s %N')
}

stop_measuring_time() {
  read e1 e2 < <(date +'%s %N')
}

show_elapsed_time() {
  s=$(bc <<< "scale=3; $((e1-s1))+($((10#$e2-10#$s2))/1000000000)")
  echo "$s seconds"
}

if [ $# -ne 1 ]; then 
  echo "Must provide a single filename as test configuration"
  exit -1
fi

filename="$1"
test_count=1
while read -r line
do
  name="$line"

  if [[ $line == "f_emboss "* ]]; then
	f_emboss=${line:9}
  fi

  if [[ $line == "f_IUPACpal "* ]]; then
	f_IUPACpal=${line:11}
  fi

  if [[ $line == "m "* ]]; then
	m=${line:1}
  fi

  if [[ $line == "M "* ]]; then
	M=${line:1}
  fi

  if [[ $line == "g "* ]]; then
	g=${line:1}
  fi

  if [[ $line == "x "* ]]; then
	x=${line:1}
  fi

  if [[ $line == "S "* ]]; then
	S=${line:1}
  fi

  if [[ $line == "" ]]; then
	echo "Running test $test_count:"
	start_measuring_time
	./IUPACpal -f $f_IUPACpal -m $m -M $M -g $g -x $x -S $S
	stop_measuring_time
	show_elapsed_time
	echo "Test $test_count complete"
	echo ""
	echo ""
	test_count=$((test_count+1))
  fi
done < "$filename"

echo "All tests completed!"

#write output to a csv file
#include reading of palindrome output analysed by python correctness and compare tools
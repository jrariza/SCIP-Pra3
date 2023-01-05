#!/bin/bash

PARTICLES=("100" "1000" "5000" "10000" "25000" "100000")
ITERATIONS=("1000" "1000" "1000" "100" "25" "10")
OUTPUT_FILE="AutoTestNBody_Output.txt"

TIMES=()
PASSED=()

filename=""
i=0

checkOutput() {
  particles=$1
  iterations=$2
  file="galaxy_""$particles""B_""$iterations""i_final.out"
  cmp=$(diff "./res/$file" "./Test/$file")
  if [ "$cmp" == "" ]; then
    echo "True"
  else
    echo "False"
  fi
}

printLine() {
  printf "%-25s %-20s %-20s %-20s %-20s %-20s %-20s %-20s\n" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" >>"$OUTPUT_FILE"
}

execSequentialMode() {
  graphics=$1

  echo "Sequential Mode"
  timeSeq=$( (time ./NBody_BarnesHut/NBody_GUI 1 "${ITERATIONS[i]}" ./NBody_BarnesHut/Test/"$filename" $graphics) 2>&1 >/dev/null |
    awk '/^real/{print $2}')
  printf "\tTiempo: %s\n" "$timeSeq"
  TIMES[1]=$timeConc
}

execConcurrentMode() {
  program=$1
  nThread=$2
  graphics=$3

  echo "Concurrent Mode ($nThread Threads)"
  timeConc=$( (time ./$program 1 -i "${ITERATIONS[i]}" -t $nThread -f "./Test/""$filename""" $graphics) 2>&1 >/dev/null |
    awk '/^real/{print $2}')
  testPassedConc=$(checkOutput "${PARTICLES[i]}" "${ITERATIONS[i]}")
  printf "\tTiempo: %s\n\tTest Passed: %s\n" "$timeConc" "$testPassedConc"
  TIMES[nThread]=$timeConc
  PASSED[nThread]=$testPassedConc
}

rm $OUTPUT_FILE
touch $OUTPUT_FILE

printLine "FILENAME" "SEQ TIME" "TIME CONC 4" "TEST PASSED 4" "TIME CONC 8" "TEST PASSED 8" "TIME CONC 16" "TEST PASSED 16"
printf "\n\nSIN Parte Gráfica:\n\n"
echo "SIN PARTE GRÁFICA" >>"$OUTPUT_FILE"

for i in {0..5}; do
  filename="galaxy_${PARTICLES[i]}B_initial.out"
  printf "\nFILE: %s\n" "$filename"

  execSequentialMode

  execConcurrentMode "NBody" 4

  execConcurrentMode "NBody" 8

  execConcurrentMode "NBody" 16

  printLine "$filename" "${TIMES[1]}" "${TIMES[4]}" "${PASSED[4]}" "${TIMES[8]}" "${PASSED[8]}" "${TIMES[16]}" "${PASSED[16]}"
done

printf "\n\nCon Parte Gráfica:\n\n"
echo "CON PARTE GRÁFICA" >>"$OUTPUT_FILE"
for i in {0..5}; do
  filename="galaxy_${PARTICLES[i]}B_initial.out"
  printf "\nFILE: %s\n" "$filename"

  execSequentialMode "-g"

  execConcurrentMode "NBody_GUI" 4 "-g"

  execConcurrentMode "NBody_GUI" 8 "-g"

  execConcurrentMode "NBody_GUI" 16 "-g"

  printLine "$filename" "${TIMES[1]}" "${TIMES[4]}" "${PASSED[4]}" "${TIMES[8]}" "${PASSED[8]}" "${TIMES[16]}" "${PASSED[16]}"
done

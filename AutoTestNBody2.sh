#!/bin/bash

PARTICLES=("100" "1000" "5000" "10000" "25000" "100000")
ITERATIONS=("1000" "1000" "1000" "100" "25" "10")
OUTPUT_FILE="AutoTestNBody_Output.txt"

checkOutput(){
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

printLine(){
  printf "%-25s %-20s %-20s %-20s %-20s %-20s %-20s %-20s\n" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" >> "$OUTPUT_FILE"
}

execConcurrentMode(){
    echo "Concurrent Mode (4 Threads)"
    timeConc4=$( (time ./NBody 1 -i "${ITERATIONS[i]}" -t 4 -f "./Test/""$filename""") 2>&1 >/dev/null |\
    awk '/^real/{print $2}')
    testPassedConc4=$(checkOutput "${PARTICLES[i]}" "${ITERATIONS[i]}")
    printf "\tTiempo: %s\n\tTest Passed: %s\n" "$timeConc4" "$testPassedConc4"
}

rm $OUTPUT_FILE
touch $OUTPUT_FILE
printLine "FILENAME" "SEQ TIME" "TIME CONC 4" "TEST PASSED 4" "TIME CONC 8" "TEST PASSED 8" "TIME CONC 16" "TEST PASSED 16"

printf "\n\nSIN Parte Gráfica:\n\n"
echo "SIN PARTE GRÁFICA">> "$OUTPUT_FILE"
for i in {1..5}; do
  filename="galaxy_${PARTICLES[i]}B_initial.out"
  printf "\nFILE: %s\n" "$filename"
#  echo "Sequential Mode"
#  timeSeq=$( (time ./NBody_BarnesHut/NBody_GUI 1 "${ITERATIONS[i]}" ./NBody_BarnesHut/Test/"$filename") 2>&1 >/dev/null |\
#  awk '/^real/{print $2}')
#  printf "\tTiempo: %s\n" "$timeSeq"

  echo "Concurrent Mode (4 Threads)"
  timeConc4=$( (time ./NBody 1 -i "${ITERATIONS[i]}" -t 4 -f "./Test/""$filename""") 2>&1 >/dev/null |\
  awk '/^real/{print $2}')
  testPassedConc4=$(checkOutput "${PARTICLES[i]}" "${ITERATIONS[i]}")
  printf "\tTiempo: %s\n\tTest Passed: %s\n" "$timeConc4" "$testPassedConc4"

  echo "Concurrent Mode (8 Threads)"
  timeConc8=$( (time ./NBody 1 -i "${ITERATIONS[i]}" -t 8 -f ./Test/"$filename") 2>&1 >/dev/null |\
  awk '/^real/{print $2}')
  testPassedConc8=$(checkOutput "${PARTICLES[i]}" "${ITERATIONS[i]}")
  printf "\tTiempo: %s\n\tTest Passed: %s\n" "$timeConc8" "$testPassedConc8"

  echo "Concurrent Mode (16 Threads)"
  timeConc16=$( (time ./NBody 1 -i "${ITERATIONS[i]}" -t 16 -f ./Test/"$filename" ) 2>&1 >/dev/null |\
  awk '/^real/{print $2}')
  testPassedConc16=$(checkOutput "${PARTICLES[i]}" "${ITERATIONS[i]}")
  printf "\tTiempo: %s\n\tTest Passed: %s\n" "$timeConc16" "$testPassedConc16"

  printLine "$filename" "$timeSeq" "$timeConc4" "$testPassedConc4" "$timeConc8" "$testPassedConc8" "$timeConc16" "$testPassedConc16"
done


printf "\n\nCon Parte Gráfica:\n\n"
echo "CON PARTE GRÁFICA">> "$OUTPUT_FILE"
for i in {5..5}; do
  filename="galaxy_${PARTICLES[i]}B_initial.out"
  printf "\nFILE: %s\n" "$filename"
  echo "Sequential Mode"

  timeSeq=$( (time ./NBody_GUI 1 "${ITERATIONS[i]}" ./NBody_BarnesHut/Test/"$filename" 1) 2>&1 >/dev/null |\
  awk '/^real/{print $2}')
    printf "\tTiempo: %s\n" "$timeSeq"

  echo "Concurrent Mode (4 Threads)"
  timeConc4=$( (time ./NBody_GUI 1 -i "${ITERATIONS[i]}" -t 4 -f ./Test/"$filename" -g) 2>&1 >/dev/null |\
  awk '/^real/{print $2}')
  testPassedConc4=$(checkOutput "${PARTICLES[i]}" "${ITERATIONS[i]}")
  printf "\tTiempo: %s\n\tTest Passed: %s\n" "$timeConc4" "$testPassedConc4"

  echo "Concurrent Mode (8 Threads)"
  timeConc8=$( (time ./NBody_GUI 1 -i "${ITERATIONS[i]}" -t 8 -f ./Test/"$filename" -g) 2>&1 >/dev/null |\
  awk '/^real/{print $2}')
  testPassedConc8=$(checkOutput "${PARTICLES[i]}" "${ITERATIONS[i]}")
  printf "\tTiempo: %s\n\tTest Passed: %s\n" "$timeConc8" "$testPassedConc8"

  echo "Concurrent Mode (16 Threads)"
  timeConc16=$( (time ./NBody_GUI 1 -i "${ITERATIONS[i]}" -t 16 -f ./Test/"$filename" -g) 2>&1 >/dev/null |\
  awk '/^real/{print $2}')
  testPassedConc16=$(checkOutput "${PARTICLES[i]}" "${ITERATIONS[i]}")
  printf "\tTiempo: %s\n\tTest Passed: %s\n" "$timeConc16" "$testPassedConc16"

  printLine "$filename" "$timeSeq" "$timeConc4" "$testPassedConc4" "$timeConc8" "$testPassedConc8" "$timeConc16" "$testPassedConc16"
done
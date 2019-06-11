#!/bin/env bash

echo "Compliling all results..."

rm energies.txt
echo "Config. \t Energy" >> energies.txt

while IFS='' read -r line || [[ -n "$line" ]]; do
  cd $line
  echo $line >> energies.txt
  grep "Total energy" $line-energy-run.out | grep -o -- "-[0-9]*\.[0-9]*" >> ../energies.txt
  cd ..
done < "submissions.txt"

#!/bin/bash

testdir=$1
echo "testing $testdir"

for exe in $(\ls $1 | grep -v ".[ch]"); do
  $testdir/$exe
  if [[ $? -ne 0 ]]; then
    break;
  fi
done

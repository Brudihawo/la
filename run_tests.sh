#!/bin/bash

testdir=$1
echo "testing $testdir"

for exe in $(\ls $1 | grep -v ".[ch]"); do
  $testdir/$exe
done

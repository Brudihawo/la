#!/bin/bash

kcov_dir="./coverage"
test_dir="./test"
exes=$(\ls $test_dir | grep -v "\.[ch]")

for exe in ${exes[@]}; do
  mkdir -p $kcov_dir/$exe
  kcov --exclude-pattern=log.c \
    --exclude-path=$test_dir $kcov_dir/$exe $test_dir/$exe 
done

mkdir $kcov_dir/full
kcov --merge $kcov_dir/full $kcov_dir/test_*


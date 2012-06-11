#!/bin/sh

cd $(dirname $0)

START=$(date +%s)
for testfile in ../testbin/*
do
    echo "Running..."
    echo $testfile
    ./$testfile
    if [ ! $? = 0 ]; then
        echo "Test failure in " $testfile
        exit -1
    fi
done

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "All tests successfully run in approximately $DIFF seconds."

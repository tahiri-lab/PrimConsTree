#!/bin/bash

tqd=$1
t1=$2
t2=$3

echo "$t1" > "temp1"
echo "$t2" > "temp2"
result=$($tqd "temp1" "temp2")
rm "temp1"
rm "temp2"

echo "$result"
#!/bin/bash

rm temp.txt > /dev/null 2>&1

$(for i in $(seq $1 $2); do python -m primconstree -f $3 -s $i -v 1 >> temp.txt; done)
cat temp.txt | sort | uniq -c

# cat temp.txt | while IFS= read -r l; do
# 	echo "$(grep -c $l temp.txt) $l" >> temp2.txt
# done 

rm temp.txt

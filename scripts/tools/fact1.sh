#!/bin/bash

output=$($1 <<EOF
$2
$3
1
EOF
)
echo "$output" | head -n 6 | tail -n 1
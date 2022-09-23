#!/bin/bash

# prime
max_cnt=100

for ((i = 0; i < $max_cnt; i++)); do
	dividend=$(($RANDOM % 97))   # 0 to 96
	divisor=$(($RANDOM % 96 + 1)) # 1 to 96
	
	echo "$dividend / $divisor = $(($dividend/$divisor))"
	
	./build/Okada_2019 "$dividend" "$divisor"
	#sleep 3
done

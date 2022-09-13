#!/bin/bash

# prime
max=17

for ((i = 0; i < $max; i++)); do
	for ((j = 1; j < $max; j++)); do
		echo "$i / $j = $(($i / $j))"
		./build/Okada_2019 $i $j
	done
done

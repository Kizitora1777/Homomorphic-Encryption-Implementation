#!/bin/bash

# prime
max=7

for ((i = 0; i < $max; i++)); do
	for ((j = 1; j < $max; j++)); do
		echo "$i / $j = $(($i / $j))"
		./Okada_2019 $i $j
	done
done

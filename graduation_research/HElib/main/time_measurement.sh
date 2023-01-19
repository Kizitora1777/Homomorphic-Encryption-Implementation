#!/bin/bash

# input size
input_size=21

# input
inputs=(-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 )

for ((i = 0; i < $input_size; ++i)) do
    echo  "x = ${inputs[i]}"
    ./build/main ${inputs[i]}
done

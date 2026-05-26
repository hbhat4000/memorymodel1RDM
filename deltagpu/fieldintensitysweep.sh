#!/bin/bash

for ((i=1;i<=10;i++)); do
    for ((j=1;j<=10;j++)); do
        for ((k=0;k<=2;k++)); do

            # delay = 100 + 50*k
            delay=$((100 + 50*k))

            # intensity: logspace(1e-4, 1e-2, 10)
            intensity=$(awk -v j="$j" 'BEGIN {
                start = 1e-4;
                end   = 1e-2;
                ratio = exp((log(end/start))/9);
                printf "%.8f", start * (ratio^(j-1));
            }')

            # freq: logspace(0.001, 0.02, 10)
            freq=$(awk -v i="$i" 'BEGIN {
                start = 1e-3;
                end   = 2e-2;
                ratio = exp((log(end/start))/9);
                printf "%.8f", start * (ratio^(i-1));
            }')

            echo $delay $intensity $freq

            ./memoryFO \
                --time 0.05,200.0 \
                --field $freq,$intensity,1 \
                --delay $delay \
                --infile ../psi4data/fci_h2_6-31gss.npz \
                --outpath ./sweep/ \
                --verbose --savetraj

        done
    done
done

for ((delay=40;delay<=760;delay+=40)); do
    ./memoryFO --time 0.01,100.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_6-31g.npz --outpath ./finescale/ --g0 0 --g1 1 --verbose
done

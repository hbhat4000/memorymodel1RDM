export CUDA_VISIBLE_DEVICES=1
for ((delay=50;delay<=950;delay+=50)); do
    ./memoryFO --time 0.005,100.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_6-31g.npz --outpath ./finescale/ --g0 0 --g1 1 --verbose
done

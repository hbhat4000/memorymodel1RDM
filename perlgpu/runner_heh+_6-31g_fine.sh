export CUDA_VISIBLE_DEVICES=2
for ((delay=950;delay<=1250;delay+=50)); do
    ./memoryFO --time 0.004,100.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_heh+_6-31g.npz --outpath ./finescale/ --g0 0 --g1 1 --verbose 
done


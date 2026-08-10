for ((delay=1;delay<=10;delay+=1)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/cis_h2_sto-3g.npz --outpath ./cisresults/sto-3g/ --g0 0 --g1 1 --verbose 
done

for ((delay=1;delay<=10;delay+=1)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/cis_h2_sto-3g.npz --outpath ./cisresults/sto-3g/ --g0 0 --g1 2 --verbose 
done

for ((delay=1;delay<=10;delay+=1)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/cis_h2_sto-3g.npz --outpath ./cisresults/sto-3g/ --g0 1 --g1 2 --verbose 
done


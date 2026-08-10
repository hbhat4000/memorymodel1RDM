for ((delay=1;delay<=10;delay+=1)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/cis_heh+_6-31g.npz --outpath ./cisresults/6-31g/ --g0 0 --g1 1 --verbose 
done

for ((delay=1;delay<=10;delay+=1)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/cis_heh+_6-31g.npz --outpath ./cisresults/6-31g/ --g0 0 --g1 2 --verbose 
done

for ((delay=1;delay<=10;delay+=1)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/cis_heh+_6-31g.npz --outpath ./cisresults/6-31g/ --g0 1 --g1 2 --verbose 
done

for ((delay=1;delay<=10;delay+=1)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/cis_heh+_6-31g.npz --outpath ./cisresults/6-31g/ --g0 0 --g1 3 --verbose 
done

for ((delay=1;delay<=10;delay+=1)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/cis_heh+_6-31g.npz --outpath ./cisresults/6-31g/ --g0 1 --g1 3 --verbose 
done

for ((delay=1;delay<=10;delay+=1)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/cis_heh+_6-31g.npz --outpath ./cisresults/6-31g/ --g0 2 --g1 3 --verbose 
done

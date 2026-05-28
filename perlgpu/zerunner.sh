for ((delay=5;delay<=955;delay+=50)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_3-21g.npz --outpath ./3-21g/ --g0 0 --g1 1 --verbose 
done

for ((delay=5;delay<=955;delay+=50)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_3-21g.npz --outpath ./3-21g/ --g0 0 --g1 2 --verbose 
done

for ((delay=5;delay<=955;delay+=50)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_3-21g.npz --outpath ./3-21g/ --g0 1 --g1 2 --verbose 
done

for ((delay=5;delay<=955;delay+=50)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_3-21g.npz --outpath ./3-21g/ --g0 0 --g1 3 --verbose 
done

for ((delay=5;delay<=955;delay+=50)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_3-21g.npz --outpath ./3-21g/ --g0 1 --g1 3 --verbose 
done

for ((delay=5;delay<=955;delay+=50)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_3-21g.npz --outpath ./3-21g/ --g0 2 --g1 3 --verbose 
done
for ((delay=5;delay<=455;delay+=50)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_heh+_3-21g.npz --outpath ./3-21g/ --g0 0 --g1 1 --verbose 
done

for ((delay=5;delay<=455;delay+=50)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_heh+_3-21g.npz --outpath ./3-21g/ --g0 0 --g1 2 --verbose 
done

for ((delay=5;delay<=455;delay+=50)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_heh+_3-21g.npz --outpath ./3-21g/ --g0 1 --g1 2 --verbose 
done

for ((delay=5;delay<=455;delay+=50)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_heh+_3-21g.npz --outpath ./3-21g/ --g0 0 --g1 3 --verbose 
done

for ((delay=5;delay<=455;delay+=50)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_heh+_3-21g.npz --outpath ./3-21g/ --g0 1 --g1 3 --verbose 
done

for ((delay=5;delay<=455;delay+=50)); do
    ./memoryFO --time 0.01,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_heh+_3-21g.npz --outpath ./3-21g/ --g0 2 --g1 3 --verbose 
done

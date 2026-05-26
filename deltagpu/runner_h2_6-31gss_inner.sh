# for ((delay=110;delay<=150;delay+=5)); do
#     ./memoryFO --time 0.05,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_6-31gss.npz --outpath ./ --g0 0 --g1 2 --verbose
# done
# for ((delay=160;delay<=200;delay+=5)); do
#     ./memoryFO --time 0.05,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_6-31gss.npz --outpath ./ --g0 0 --g1 2 --verbose
# done

# MEMORY ERRORS?
# for ((delay=110;delay<=150;delay+=5)); do
#     echo $delay 1 2
#     ./memoryFO --time 0.05,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_6-31gss.npz --outpath ./ --g0 1 --g1 2 --verbose
# done

# MEMORY ERRORS?
# for ((delay=160;delay<=200;delay+=5)); do
#     echo $delay 1 2
#     ./memoryFO --time 0.05,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_6-31gss.npz --outpath ./ --g0 1 --g1 2 --verbose
# done

# MEMORY ERRORS?
# for ((delay=110;delay<=150;delay+=5)); do
#     echo $delay 0 3
#     ./memoryFO --time 0.05,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_6-31gss.npz --outpath ./ --g0 0 --g1 3 --verbose
# done

# MEMORY ERRORS?
# for ((delay=160;delay<=200;delay+=5)); do
#     echo $delay 0 3
#     ./memoryFO --time 0.05,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_6-31gss.npz --outpath ./ --g0 0 --g1 3 --verbose
# done

for ((delay=110;delay<=150;delay+=5)); do
    echo $delay 1 3
    ./memoryFO --time 0.05,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_6-31gss.npz --outpath ./ --g0 1 --g1 3 --verbose
done

for ((delay=160;delay<=200;delay+=5)); do
    echo $delay 1 3
    ./memoryFO --time 0.05,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_6-31gss.npz --outpath ./ --g0 1 --g1 3 --verbose
done

for ((delay=110;delay<=150;delay+=5)); do
    echo $delay 2 3
    ./memoryFO --time 0.05,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_6-31gss.npz --outpath ./ --g0 2 --g1 3 --verbose
done

for ((delay=160;delay<=200;delay+=5)); do
    echo $delay 2 3
    ./memoryFO --time 0.05,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_h2_6-31gss.npz --outpath ./ --g0 2 --g1 3 --verbose
done

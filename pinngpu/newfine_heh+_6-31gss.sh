export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/u/hbhat/lib:/u/hbhat/lib64:/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/math_libs/12.8/lib64
export CUDA_VISIBLE_DEVICES="1"

for ((delay=10;delay<=1110;delay+=10)); do
    # Construct the expected output filename
    target_file="./6-31gss/fci_heh+_6-31gss_0.025000_${delay}_0.127827_0.100000_5.txt"

    # Check if the file exists; if it does, skip the rest of this iteration
    if [ -f "$target_file" ]; then
        echo "Detected existing file: $target_file"
        continue
    fi
    ./memoryFO --time 0.025,200.0 --field -1,0.1,5 --delay $delay --infile ../psi4data/fci_heh+_6-31gss.npz --outpath ./6-31gss/ --g0 0 --g1 1 --verbose
done

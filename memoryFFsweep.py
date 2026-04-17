import numpy as np
import os
import subprocess
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(prog='memoryFFsweep',
                                 description='Record the MAE (mean absolute error) of the 1RDM field-free memory model as we sweep through values of the delay')

parser.add_argument('--dt', type=float, required=True, help='time step')
parser.add_argument('--mindelay', type=int, required=True, help='minimum delay')
parser.add_argument('--maxdelay', type=int, required=True, help='maximum delay')
parser.add_argument('--delaystep', type=int, required=True, help='delay step')
parser.add_argument('--infile', required=True, help='input file')

# actually parse command-line arguments
args = parser.parse_args()

# time step
dt = args.dt

# delay range
numdelays = ((args.maxdelay - args.mindelay) // args.delaystep) + 1
delayrange = args.mindelay + np.arange(numdelays,dtype=np.int16)*args.delaystep
maes = np.zeros(numdelays)

# input file
infile = args.infile

# Run the command and capture output
j = 0
for thisdelay in delayrange:
    result = subprocess.run(['./memoryFF', '--dt', str(dt), '--delay', str(thisdelay), '--infile', str(infile)], capture_output=True, text=True)
    maes[j] = float(result.stdout)
    print("delay = " + str(thisdelay) + "; mae = " + str(maes[j]))
    j += 1    

p = Path(infile)
stem = p.stem
outfile = 'maes_' + stem + '.npz'

# Save output
np.savez(outfile, maes=maes, delayrange=delayrange, dt=dt)



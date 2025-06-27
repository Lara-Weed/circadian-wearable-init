#!/bin/bash

# List your subject numbers here:
SUBS=( $(seq 1 1 88) )

for subs in "${SUBS[@]}"; do
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=CAs${subs}
#SBATCH --output=CAs${subs}.%j.out
#SBATCH --error=CAs${subs}.%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --mail-type=ALL

# Load MATLAB
module load matlab

# Respect the requested core count in MATLAB's parallel pool
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK

echo "Submitting job for Sub = ${subs}"

matlab -batch "parpool('local', \$SLURM_CPUS_PER_TASK); runLightModel_Wearable_ODE15_AMP(${subs}); exit"
EOF

done


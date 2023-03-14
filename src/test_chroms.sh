#!/bin/bash 

for i in {1..22} X
do
for cutoff in -0.1 -1.0 -2.0 -5.0
do
for window_size in 50 500 1000 2000 4000 8000 12000
do
job_script="test_chr"$i"_"${window_size}"_"${cutoff}
echo "#!/bin/bash" >> ${job_script}
echo "#SBATCH -c 1" >> ${job_script}                               # Request one core
echo "#SBATCH -N 1" >> ${job_script}                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
echo "#SBATCH --mem=100G" >> ${job_script}                          # Memory total in MB (for all cores)
echo "#SBATCH -o slurm_%j.out" >> ${job_script}                 # File to which STDOUT will be written, including job ID
echo "#SBATCH -e slurm_%j.err" >> ${job_script}                 # File to which STDERR will be written, including job ID

echo	"./linker pop -g /czlab/alumni/tourdot/cpp_linker_v3.6_HCC1954/output/graph_variant_oct23_BL1954_all_hic_chr"${i}".dat -v /czlab/gbrunette/Eagle_v2.4.1/HCC1954BL_21-6-4/HCC1954BL_chr"${i}"_21-6-4.vcf.gz -e "${cutoff}" -w "${window_size}" -c chr"$i" -n test_"${window_size}"_"${cutoff}" >> ${job_script}
done
done
done

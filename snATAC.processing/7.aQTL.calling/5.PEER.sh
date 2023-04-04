for ct in Ast Ex In Microglia OPC Oligo  PerEndo
do
echo "#!/bin/bash
#SBATCH -N 1                                # Request one node
#SBATCH -n 6                               # Request n cpu (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 1-12:00                          # Runtime in D-HH:MM format
#SBATCH -o $ct.peer.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e $ct.peer.err                 # File to which STDERR will be written, including job ID
#SBATCH -p kellis
#SBATCH --mem=60000                         ## member used
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=xiongxushen1992@gmail.com    # Email to which notifications will be sent

module load r/3.5.1
Rscript run_PEER.R $ct.TMM.voom.normed.mtx.sort.bed.gz $ct.TMM.voom.normed 50" > $ct.runPEER.sh
sbatch $ct.runPEER.sh
done


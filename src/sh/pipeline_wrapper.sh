############################
# Wrapper to run Snakemake #
############################

# The following options are included in every command because they are generally useful
generalOptions="--rerun-incomplete --latency-wait 60 --keep-going"

## explanation: 
## Re-run all jobs the output of which is recognized as incomplete, wait this time before considering a job does not have output, keep going with other jobs if one fails

# Options helpful for debugging
debuggingOptions="--printshellcmds --verbose --reason"

## explanation:
## print out the shell commands that will be executed, the debugging output and the reason for each rule to be executed

## It is required that genome indexing has already been performed on login.cluster.embl.de using:
## sbatch /g/boulard/Projects/BS_EMSeq/src/sh/indexing_genome.sh
## cp -r Bisulfite_Genome/ /g/boulard/Projects/BS_EMSeq/data/genome/

# Going to the main folder of the experiment folder, which will be the parent of all folders in the pipeline and in the Snakefile
cd /g/boulard/Projects/BS_EMSeq/

snakemake \
  -s snake-make/Snakefile \
  --configfile config/config.yaml \
  --use-conda \
  $generalOptions\
  $debuggingOptions\
  --use-singularity --singularity-args "--bind /g/boulard/Projects/BS_EMSeq/" \
  -j 20 \
  --cluster-config config/cluster.SLURM.json \
  --cluster "sbatch -p {cluster.queue} -J {cluster.name} --cpus-per-task {cluster.nCPUs} --mem {cluster.memory} \
                   --time {cluster.maxTime} -o {cluster.output} -e {cluster.error}" \
  --local-cores 1


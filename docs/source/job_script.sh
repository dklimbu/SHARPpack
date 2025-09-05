#!/bin/bash
#============================================================
#     JOB SCRIPT TO RUN SHARP PACK v2 ON HPC SERIAL/PARALLEL
#     WHILE RUNNING PARALLEL, BASED ON n-CPU, IT CREATES
#     n-DIRECTORIES AND RUN JOBS PARALLELLY, AND FINALLY 
#     COMPUTES AVERAGE OF THE DISIRED POPULATATION (NEED TO
#     SPECIFY FOR INPUT, SEE EXAMPLE IN, WRITE_AVERAGE_INPUT)
#
#     authors    - D.K. Limbu & F.A. Shakib     
#
#     Method Development and Materials Simulation Laboratory
#     New Jersey Institute of Technology
#
#     USAGE:: bash job_script_hpc.sh
#           - creates and run job(s)
#
#============================================================

# Path to your executable
#============================================================
root2bin=~/Softwares/SHARP_pack2/bin
exe=${root2bin}/sharp.x

#============================================================
maxcore=1
ncore=$(awk '/ncore/{print $2}' param.in)
maxcore=$((ncore > maxcore ? ncore : maxcore))
echo 'ncore(s):' $maxcore

# ==== CONFIGURATION ====
#============================================================
root_dir="run"            # Base directory name
ncore=$maxcore            # Number of parallel jobs/directories
param_file="param.in"     # Parameter file to copy to each directory

# Determine padding width based on ncore
if [ "$ncore" -lt 10 ]; then
   pad=1
elif [ "$ncore" -lt 100 ]; then
   pad=2
else
   pad=3
fi

# ==== FUNCTION TO CREATE AVERAGE INPUT ====
write_average_input() {
  #fname :: File name of specific population data for averaging
  #df    :: No. of data column(s)(except first column)?

  fname="pop_diabat3.out"   
  df=2

  echo "Calculating average results from $ncore parallel jobs."

  rm -f input

  echo ${root_dir} >> input
  echo ${fname} >> input
  echo ${ncore} >> input
  echo ${df} >> input
  echo ${pad} >> input

}

# ==== FUNCTION TO CREATE PARALLEL DIRECTORIES ====
paralleldir() {
  # Check if the first directory already exists
  first_dir=$(printf "${root_dir}%0${pad}g" 1)
  if [ -d "$first_dir" ]; then
    echo '**********************************'
    echo "$first_dir exists. Remove it first!"
    echo '**********************************'
    exit 1
  fi

  # Create directories and copy files
  for i in $(seq -f "%0${pad}g" 1 $ncore); do
    dir="${root_dir}${i}"
    mkdir "$dir"
    cp param.in "$dir/param.in"
    echo "$i" > "$dir/icpu.in"
  done

  echo "Created $ncore directories: ${root_dir}$(printf "%0${pad}d" 1) to ${root_dir}${ncore}"
}

# ==== FUNCTION TO WRITE SLURM ARRAY SCRIPT ====
write_array_job() {

  cat <<EOF > submit_array.sh
#!/bin/bash
#SBATCH --job-name=array_run
#SBATCH --array=1-${ncore}
#SBATCH --output=${root_dir}%0${pad}a/myoutput_%A_%a.out
#SBATCH --error=${root_dir}%0${pad}a/jobError_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p RM-shared
#SBATCH --mem-per-cpu=2000M
#SBATCH --mail-type=end,fail
#SBATCH --time=01:00:00

cd ${root_dir}\$(printf "%0${pad}d" \${SLURM_ARRAY_TASK_ID})

 $exe

EOF

  chmod +x submit_array.sh
}

# ==== FUNCTION TO WRITE SERIAL SLURM  SCRIPT ====
write_serial_job() {

  cat <<EOF > submit_serial.sh
#!/bin/bash
#SBATCH --job-name=serial_run
#SBATCH --output=myoutput%j.out
#SBATCH --error=jobError.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p RM-shared
#SBATCH --mem-per-cpu=2000M
#SBATCH --mail-type=end,fail
#SBATCH --time=10:00:00

 $exe

EOF

  chmod +x submit_serial.sh
}

# ==== FUNCTION TO WRITE ANALYSIS JOB SCRIPT ====
write_average_job() {
  cat <<EOF > average_job.sh
#!/bin/bash
#SBATCH --job-name=average
#SBATCH --output=myoutput%j.out
#SBATCH --error=jobError.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p RM-shared
#SBATCH --mem-per-cpu=2000M
#SBATCH --mail-type=end,fail
#SBATCH --time=00:10:00

  ${root2bin}/average.x

EOF

  chmod +x average_job.sh
}

#============================================================
## RUN SERIAL JOB ##
#============================================================
if [ $ncore -eq 1 ] ; then

  echo "Serial job running in a single core!!"

  write_serial_job

  sbatch submit_serial.sh

#============================================================
## RUN PARALLEL JOB(S) ON HPC ##
#============================================================
elif [ $ncore -gt 1 ] ; then
  #create parallel jobs

  echo "Parallel job running on " $ncore "core(s)."

  # ==== MAIN ====
  paralleldir
  write_array_job
  write_average_input
  write_average_job

  # Submit array job
  array_jobid=$(sbatch --parsable submit_array.sh)

  echo "Submitting post-processing job dependent on array job ${array_jobid}..."

  # Submit average job to run only after array job completes
  sbatch --dependency=afterok:$array_jobid average_job.sh

  echo "All jobs submitted successfully!"

fi

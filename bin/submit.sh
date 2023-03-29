#!/bin/sh
# JOB SCRIPT TO RUN SHARP PACK
#
# DIL LIMBU
# 
# jobId=1 -> running parallel jobs
# jobId=2 -> averaging parallel job results

jobId=1


root2bin=~/Desktop/MyDocs/SHARPpack/bin
exe=${root2bin}/sharp.x

#====================================================
maxcpu=1
ncpu=$(awk '/ncpu/{print $2}' param.in)

maxcpu=$((ncpu > maxcpu ? ncpu : maxcpu))

ncpu=$maxcpu

echo 'cpu:' $ncpu

root_dir=run

paralleldir(){
  if [ $ncpu -lt 10 ];
  then
    DIR=${root_dir}1
  else
    DIR=${root_dir}01
  fi

  if [ ! -d ${DIR} ]
  then
    if [ $ncpu -lt 10 ]
    then
      mkdir $(seq -f "${root_dir}%01g" 1 $((ncpu)))
    else
      mkdir $(seq -f "${root_dir}%02g" 1 $((ncpu)))
    fi
    echo ${root_dir} created!!
  else
    echo '**********************************'
    echo  ${DIR} exists, remove it first!! 
    echo 'OR Calculate average with jobID=2'
    echo '**********************************'
    echo
    exit 
  fi

  for folder in ${root_dir}*
  do
    cd $folder
      
      cp ../param.in param.in

      echo $folder | grep -Po "${root_dir}\K\d+" > icpu.in
  
      #nohup ../job.sh &
  
      ${exe} & 

    cd ..

  done
}

if [ $ncpu == 1 ];then

  ${exe} &

  echo "Serial job running in a single cpu!!"

elif [ $ncpu -gt 1 ] && [ $jobId == 1 ] ;then
  #create parallel jobs
  paralleldir

  echo "Parallel job running on " $ncpu "cpus."

elif [ $ncpu -gt 1 ] && [ $jobId == 2 ] ;then
  if [ $ncpu -lt 10 ];
  then
    DIR=${root_dir}1
    iformat=1
  else
    DIR=${root_dir}01
    iformat=2
  fi

  if [ ! -d ${DIR} ]
  then
    echo '**********************************'
    echo 'PARALLEL JOB RUN DOES NOT EXIST!!'
    echo 'USE jobID=1 to run parallel job(s)'
    echo "at FIRST with ncpu > 1 in param.in"
    echo "**********************************"
    echo "  " 
    exit
  fi

  rm -f input

  echo "Calculating average results from " $ncpu "parallel jobs."

  echo "Enter file name:"
  read fname
  echo "No. of data column(s)(except first column)?"
  read df

  echo ${root_dir} >> input
  echo ${fname} >> input
  echo ${ncpu} >> input
  echo ${df} >> input
  echo ${iformat} >> input

  ${root2bin}/average.x

  rm -f input
else

  echo "Error in jobId!!, USE"
  echo "jobId=1 for running parallel jobs"
  echo "jobId=2 for averaging parallel job results"
  echo 
fi


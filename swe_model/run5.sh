
#set number of nodes and processes
PBS -l nodes=1 ppn=20
#set max wall clock time
PBS -l walltime=900:00:00
#set name of the job
PBS -N run5
#mail alert at the start end and abortion of execution
PBS -m bea
#send mail to this address
PBS -M hkropp@colgate.edu
#use submission environment
PBS -V
#start job from the directory it was submitted
cd $PBS_0_WORKDIR
#run the script
/shared/R-3.4.3/bin/R Rscript /home/hkropp/github/boreal_lw/swe_model/swe_depletion_model_script_r5.r

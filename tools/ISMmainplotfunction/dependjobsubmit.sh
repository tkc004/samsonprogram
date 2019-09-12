declare -a arr=("element1" "element2" "element3")
RES=$(qsub job.sh "${arr[0]}" )
sbatch --dependency=afterany:${RES##* } job.sh "${arr[1]}"
sbatch --dependency=afterany:${RES##* } job.sh "${arr[2]}"
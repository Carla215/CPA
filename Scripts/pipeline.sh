#!/bin/bash

#load modules
module load anaconda3/personal
module load ldsc
#ldsc13 is a specific environment with the right LDSC dependencies (runs impossibly slowly without this)
source activate ldsc13
echo "LDSC loaded"

########
# LDSC #
########

#get command-line arguments for traits and cohort sizes, using flags -t and -n
while getopts ":t:n:" opt; do
  case $opt in
    t) set -f
	trait="$OPTARG"
      ;;
    n) set -f
	IFS=','
	number="$OPTARG"
      ;;
    \?) echo "Usage: cmd [-t] [-n]"
      ;;
  esac
done
shift $((OPTIND-1))



sumstat_list=()
#split command-line options into separate variables, add to list
for var in "${trait[@]}";
do
	IFS=',' read -a traits <<< "$trait"
	IFS=',' read -a numbers <<< "$number"
done


# LDSC: munge summary statistics
for var in ${!numbers[*]};
do
	munge_sumstats.py \
	--sumstats "${traits[$var]}" \
	--N "${numbers[$var]}" \
	--out "${traits[var]}" \
	--merge-alleles w_hm3.noMHC.snplist
done

#add .sumstats.gz to end of trait names
for var in "${!traits[*]}";
do
	sumstat_list+=($var.sumstats.gz)
done

#join by commas
join() {
    local retname=$1 sep=$2 ret=$3
    shift 3 || shift $(($#))
    printf -v "$retname" "%s" "$ret${@/#/$sep}"
}

join sumstats "," "${sumstat_list[@]}"


#LDSC: ldsc (one-vs-all; must update to all-vs-all)
ldsc.py \
--rg "$sumstats" \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out ldsc_"$trait"

#deactivate ldsc environment
conda deactivate

#activate environment with R dependencies
source activate Rker

# R script containing GenomicSEM & ASSET commands
Rscript pipeline.R $trait $number






#!/bin/bash

# getopts
OPTS=`getopt -o CELAP:M:Nh --long conserved,errors,lengths,abundences,project:,modifier:,picks,help -n 'parse-options' -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

conserved="false"
errors="false"
lengths="false"
project="default"
modifier="simulator"
picks="false"

while true; do
  case "$1" in
  -C | --conserved ) conserved="true"; shift ;;
  -E | --errors ) errors="true"; shift ;;
  -L | --lengths ) lengths="true"; shift ;;
  -A | --abundences ) abundences="true"; shift ;;
  -P | --project ) project="$2"; shift 2;;
  -M | --seq-mod ) modifier="$2"; shift 2;;
  -N | --picks ) picks="true"; shift;;
  -h | --help ) help=true; shift;;
  -- ) shift; break ;;
  * ) break ;;
  esac
done
shift $((OPTIND -1))

if [ $help ]; then
 echo "-C or --conserved - add conserved sequences to all OTUs"
 echo "-E or --errors - add errors to all OTUs"
 echo "-L or --lengths - add length variance to OTUs"
 echo "-A or --abundances - add High, Middling, Low abundances to OTUs"
 echo "-M or --seq-mod - add a sequence modifier"
 echo "-P or --project - add a project name so that there are not collisions"
 exit 1
fi

region="CGTCACACTTCATGATGGAATTGA"

min_max_lengths=" -L 500 -l 500"
if [ $lengths == "true" ]; then min_max_lengths=" -L 500 -l 350"; fi

min_max_errors=" -D 0 -d 0"
if [ $errors == "true" ]; then min_max_errors=" -D 10 -d 1"; fi

min_max_picks=" -W 136 -w 136"
if [ $picks == "true" ]; then min_max_picks=" -W 1360 -w 14"; fi

min_max_picks_abun=" -W 45 -w 45"
if [ $picks == "true" ]; then min_max_picks_abun=" -W 453 -w 5"; fi

num_fasta=68
num_OTU=68

if [ $abundences == "true" ]; then
  echo -e "abundances\n"
  ./Barcode_Simulator -O -o $(($num_OTU/6)) -s 10 $min_max_errors $min_max_lengths -p $project -S "$project-high-otus"
  ./Barcode_Simulator -O -o $(($num_OTU/6*2)) -s 10 $min_max_errors $min_max_lengths -p $project -S "$project-middling-otus"
  ./Barcode_Simulator -O -o $(($num_OTU/6*3)) -s 10 $min_max_errors $min_max_lengths -p $project -S "$project-low-otus"
else
  echo -e "no abundances\n"
  ./Barcode_Simulator -O -o $num_OTU -s 10 $min_max_errors $min_max_lengths -p $project -S "$project-otus"
fi

if [ $conserved == "true" ]; then
 if [ $abundences == "true" ]; then
  for abun in "high middling low"; do
   for line in $( cat "$project-$abun-otus" ); do
   echo -e "$region$line"
   done >> "$project-otus-2"
   mv "$project-otus-2" "$project-$abun-otus"
  done
 else
  for line in $( cat "$project-otus" ); do
   echo -e "$region$line"
  done >> "$project-otus-2"
  mv "$project-otus-2" "$project-otus"
 fi
fi

if [ $abundences == "true" ]; then
 ./Barcode_Simulator -p $project -f $num_fasta -r "$project-high-otus" $min_max_picks_abun -a "-$modifier-high"
 ./Barcode_Simulator -p $project -f $num_fasta -r "$project-middling-otus" $min_max_picks_abun -a "-$modifier-middling"
 ./Barcode_Simulator -p $project -f $num_fasta -r "$project-low-otus" $min_max_picks_abun -a "-$modifier-low"
else
 ./Barcode_Simulator -p $project -f $num_fasta -r "$project-otus" $min_max_picks -a "-$modifier"
fi

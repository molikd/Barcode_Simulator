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
modifier="sequence"
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

#./Run_Simulation.sh
DIR=$( pwd )
CORES=16

options=$( 
if [[ $conserved == "true" ]]; then echo -n " -C"; fi
if [[ $errors == "true" ]]; then echo -n " -E"; fi
if [[ $lengths == "true" ]]; then echo -n " -L"; fi 
if [[ $abundences == "true" ]]; then echo -n " -A"; fi
if [[ $picks == "true" ]]; then echo -n " -N"; fi 
)

echo "Selected: $options"

./Run_Simulation.sh $options --project $project --modifier $modifier

FASTAS=$( ls | grep -E "^$project-[0-9]+\.fasta" )

for file in $FASTAS; do cat $file >> $project-all.fasta; done

#pick_de_novo_otus.py -a -O $CORES -i $DIR/default_all.fasta -o $DIR/def-$project -p $DIR/qiime_params.txt

pick_open_reference_otus.py -a --min_otu_size 1 -n denovo --suppress_align_and_tree --force -o $DIR/${project}-def -O $CORES -r $DIR/${project}-ref.fasta -i "$( echo $FASTAS | sed -e 's/ /,/g' )" 

biom convert -i $DIR/${project}-def/otu_table_mc1_w_tax.biom -o $DIR/results/${project}_otu_table.txt --to-tsv --table-type="OTU table"

for file in $FASTAS; do 
 mash sketch -p $CORES $file
done

for file in $(  ls | grep "$project-[0-9]*\.fasta.msh" ); do 
 mash dist -p $CORES $file $( ls | grep "$project-[0-9]*\.fasta.msh" )
done | awk '{ if ($1 < $2) { print $1" "$2" "$3 } else { print $2" "$1" "$3 } }' | sort | uniq | grep -v " 0$" >> $DIR/results/${project}_mash_dists.txt

rm ${project}-ref.fasta
rm ${project}-otus
rm ${project}-all.fasta
rm -rf ${project}-def
for file in $FASTAS; do rm $file; done
for file in $( ls | grep -E "^$project-[0-9]+\.fasta.msh" ); do rm $file; done

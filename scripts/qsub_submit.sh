#!/bin/bash

#echo "--project"

combine() {
    local limit=$[ 1 << $# ]
    local args=("$@")
    for ((value = 1; value < limit; value++)); do
        local parts=()
        for ((i = 0; i < $#; i++)); do
            [ $[(1<<i) & value] -ne 0 ] && parts[${#parts[@]}]="${args[i]}"
        done
        echo -n "${parts[@]}"
		echo
    done
}

IFS=$'\n' 
#for option in $( combine --conserved --errors --lengths --abundences --picks ); do 
for option in $( combine -C -E -L -A -N); do
 opt=$( echo $option | sed -e 's/-//g' -e 's/ //g' )
 cp qsub_simulator_pipeline_template.sh qsub_simulator_pipeline_template_temp.sh
 sed -i "s/simulator_options/$option/g" qsub_simulator_pipeline_template_temp.sh
 sed -i "s/project_name/$opt/g" qsub_simulator_pipeline_template_temp.sh
 qsub qsub_simulator_pipeline_template_temp.sh
 rm qsub_simulator_pipeline_template_temp.sh
done

cp qsub_simulator_pipeline_template.sh qsub_simulator_pipeline_template_temp.sh
sed -i "s/simulator_options//g" qsub_simulator_pipeline_template_temp.sh
sed -i "s/project_name/NULL/g" qsub_simulator_pipeline_template_temp.sh
qsub qsub_simulator_pipeline_template_temp.sh
rm qsub_simulator_pipeline_template_temp.sh

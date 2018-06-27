#!/bin/bash

echo "File Abundance Conserved Errors Lengths Picks Type" >> meta.data.file.txt
for file in $( ls | grep "_run"); do echo $file $(
Type=""
if [[ $file = *"A"* ]]; then echo "1"; Type="${Type}A"; else echo "0"; fi
if [[ $file = *"C"* ]]; then echo "1"; Type="${Type}C"; else echo "0"; fi
if [[ $file = *"E"* ]]; then echo "1"; Type="${Type}E"; else echo "0"; fi
if [[ $file = *"L"* && $file != *"U"* ]]; then echo "1"; Type="${Type}L"; else echo "0"; fi
if [[ $file = *"N"* && $file != *"U"* ]]; then echo "1"; Type="${Type}N"; else echo "0"; fi
if [[ $Type ]]; then echo "$Type"; else echo "O"; fi
)
done >> meta.data.file.txt

for file in $( ls | grep "_run_otu_" ); do cat $file | sed -e 's/# Constructed from biom file//g' -e 's/#OTU ID/OTU.id/g' >> ${file}.1; done
for file in $( ls | grep "_run_otu_" | grep '.1' ); do mv $file $( echo $file | sed 's/\.1//g'); done

for file in $( ls | grep "_run_mash_"); do cat $file | sed -e 's/_run-/FASTA\./g' -e 's/\.fasta//g' -e 's/[A-Z]*_[0-9]*//g' >> ${file}.1; done
for file in $( ls | grep "_run_mash_" | grep '.1'); do mv $file $( echo $file | sed 's/\.1//g'); done

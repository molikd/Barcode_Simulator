#!/bin/bash
#$ -M dmolik@nd.edu
#$ -m abe
#$ -pe smp 8
#$ -N qiime_script

export PATH=/afs/crc.nd.edu/user/d/dmolik/bin:/afs/crc.nd.edu/user/d/dmolik/local/lib/perl5/bin:/afs/crc.nd.edu/user/d/dmolik/local/bin:/afs/crc.nd.edu/user/d/dmolik/local/opt/sratoolkit.2.8.2-1-centos_linux64/bin:/afs/crc.nd.edu/user/d/dmolik/local/opt/mothur:$PATH
export PATH=/afs/crc.nd.edu/user/d/dmolik/local/opt/ncbi-blast-2.6.0+/bin:/afs/crc.nd.edu/user/d/dmolik/local/opt/trout-0.9/bin:$PATH
export PYTHONPATH=/afs/crc.nd.edu/user/d/dmolik/local/lib/python2.7/site-packages
export PERL5LIB=/afs/crc.nd.edu/user/d/dmolik/local/lib/perl5/lib/perl5:/afs/crc.nd.edu/user/d/dmolik/local/lib/perl5/lib:/afs/crc.nd.edu/user/d/dmolik/local/lib/perl5/site_perl
export MANPATH=/afs/crc.nd.edu/user/d/dmolik/local/lib/perl5/man

module load bio/qiime

./Simulator_pipeline.sh

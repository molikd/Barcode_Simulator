# Barcode_Simulator
For simulating barcode data

In preparation for the K-mer/OTU comparison paper I wrote a quick and dirty OTU simulator, which I wanted for two reasons, a fine control of the kind of data that the simulator is producing, and as a learning experience to really think about the dynamics of OTU data. It can randomly vary the length, sequence to sequence differences, and number of sequences per file generated.

I wrote it in bash (of all things) and I'll attach it here.

if you want to run it on a MAC, I'd recommend installing its dependencies with brew:

brew install coreutils

brew install gnu-getopt

The code runs natively in linux.

Barcode_Simulator - the barcode simulator 
README.md - Readme file

Other Scripts:
  Barcode_Simulator_Post.R - multithreaded post analysis
  Barcode_Simulator_Post_Single.R - singlethreaded post analysis
  Run_Simulation.sh - Runs a single simulation
  Simulator_pipeline.sh - runs the post simulation analysis 
  generate_table.sh - generates a metadata file of simulations 
  qsub_simulator_pipeline_template.sh - templeate for submitting jobs
  qsub_submit.sh - submits SGE/UGE simulation job

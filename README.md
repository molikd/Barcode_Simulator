# Barcode_Simulator
For simulating barcode data

In preparation for the K-mer/MASH/NMF comparison paper I wrote a quick and dirty OTU simulator, which I wanted for two reasons, a fine control of the kind of data that the simulator is producing, and as a learning experience to really think about the dynamics of OTU data. It can randomly vary the length, sequence to sequence differences, and number of sequences per file generated.

I wrote it in bash (of all things) and I'll attach it here.

if you want to run it on a MAC, I'd recommend installing its dependencies with brew:

brew install coreutils

brew install gnu-getopt

The code runs natively in linux.

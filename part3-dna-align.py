# Part 3
# In the Data section of the course website you will find a .zip file titled DNA_Genome.zip with 15 genomes. 
#  Choose 1 of the genomes and use the dynamic programming procedure 
#  and code modules for sequence alignment that was presented in Chapter 25.6 
#  to perform sequence alignment with 4 of the other sequences. 
# Report the results of the sequence alignments you chose. 
# This is an individual assignment and should be your own work. 
import dynaprog as dyn
import blosum as b
# get 1 genome
gen1 = open('DNA_Genome/960523 01.seq').read()
#  with 4 of the other sequences. 
gen2=open('DNA_Genome/960523 02.seq')
gen3=open('DNA_Genome/960523 03.seq')
gen4=open('DNA_Genome/960523 04.seq')
gen5=open('DNA_Genome/960523 05.seq')
other_genes=[gen2,gen3,gen4,gen5]
# perform sequence alignment w/ FastNW
for gene in other_genes:
    subvals = dyn.subvals(b.BLOSUM50, b.PBET, gen1, gene)
    scormat, arrow = dyn.FastNW(subvals, gen1, gene)
    
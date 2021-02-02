from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc
import matplotlib as mpl
from matplotlib import cm
from math import log10
import math
import matplotlib.patches as patches
from matplotlib.lines import Line2D

#sequence length is 2766  bp
#molecule type is ss-DNA (single stranded DNA)
#circular viral sequences (VRL)
#features: locations/ quantifiers

# CDS is the coding sequence i.e. region of nucleotides that corresponds 
# with the sequence of amino acids in a protein

#translation - the amino acid translation - corresponds to the nucleotide coding sequence (cds)

#complement - indicates that the feature islocated on the complementary strand

#By knowing a DNA sequence, it is possible to determine the sequence of amino acids in the polypeptide chain for which the DNA codes

#first the file needs to be parsed
#get dna sequence
#seperate gene features using base span 

# gene v2 has base span 139..480
# gene v1 has base span 299..1075
# gene c3 has base span on compliment strand 1072..1476
# gene c2 has base span on compliment strand 1217..1624
# gene c1 has base span on compliment strand 1533..2612

with open("Genome.gb", "r") as handle:
     for seq_record in SeqIO.parse(handle, "genbank"):
          gen_sequence = seq_record.seq #get the genome sequence
          for feature in seq_record.features:
               if feature.type == "CDS":
                    base_span = feature.location
                    sequence = feature.location.extract(seq_record).seq
     

#using physical mapping which analyzes the physical
#distance between the known DNA sequences by calculating
#the number of base pairs between them

# gene v1 and v2 overlap
# v1 and c3 overlap
# gene c1, c2 &  c2, c3 overlap 
arr = list(gen_sequence) #put the gen sequence into an array of characters
base_names = []
base_colors = []
base_size=[]
for x in range(0,len(arr)): #go through bases and assign colors and increment size for each base
     if arr[x] == 'A':
          base_colors.append('Blue')
          base_size.append(1)
     if arr[x] == 'C':
          base_colors.append('Red') 
          base_size.append(1)  
     if arr[x] == 'G':
          base_colors.append('Green')
          base_size.append(1)
     if arr[x] == 'T':
          base_colors.append('Yellow')
          base_size.append(1)

fig, ax = plt.subplots()
names = '', 'v2', '', 'c3','','c1','' 
overlap_names = '','v1','','c2', ''
size=[139,341,592,404,57,1079,154]
my_circle=plt.Circle( (0,0), 0.8, color='white')
width = 0.15
bp_ticks_names = ['','','','','','','','','','','','300','','','','','','','','','','','','600','','','','','','','','','','','','900','','','','','','','','','','','', '1200','','','','','','','','','','', '', '1500', '','','','','','','','','','','','1800', '','','','','','','','','','','', '2100','','','','','','','','','','', '', '2400','','','','','','','','','','', '', '2700']

outside, _ = ax.pie(size,  colors=['white','darkblue','white','cornflowerblue','white','blue','white'], radius=1)
inside, _ = ax.pie([299,776,142, 407,1142], colors=['white','slateblue','white', 'blueviolet', 'white'], radius=1-width)
#ticks to use as a guide for bp where genes are located
bp_ticks, _ = ax.pie([50,3,50,3,50,3,50,3,50,3,50,10,50,3,50,3,50,3,50,3,50,3,50,10,50,3,50,3,50,3,50,3,50,3,50,10,50,3,50,3,50,3,50,3,50,3,50,10,50,3,50,3,50,3,50,3,50,3,50,10,50,3,50,3,50,3,50,3,50,3,50,10,50,3,50,3,50,3,50,3,50,3,50,10,50,3,50,3,50,3,50,3,50,3,50,10,50,3,50,3,50,3,50,3,50,3,50,10], labels=bp_ticks_names, colors=['white','black', 'white','black','white', 'black', 'white', 'black', 'white', 'black', 'white', 'black', 'white', 'black', 'white', 'black', 'white', 'black', 'white', 'black'], radius=1.2, labeldistance=1.03, textprops={'fontsize': 8})
bases, _ = ax.pie(base_size, colors=base_colors, radius=1.58)
p=plt.gcf() #get current figure

style = "Simple, tail_width=0.5, head_width=4, head_length=8"
kw = dict(arrowstyle=style, color="k")

v1 = patches.FancyArrowPatch((.5, 0.4), (-.5, 0.4),
                             connectionstyle="arc3,rad=.39", **kw)
v2 = patches.FancyArrowPatch((.75, 0.24), (.6, 0.5),
                             connectionstyle="arc3,rad=.2", **kw)
c3 = patches.FancyArrowPatch((-.74, 0.29), (-.62, 0.52),
                             connectionstyle="arc3,rad=-.2", **kw)
c2 = patches.FancyArrowPatch((-.54, -0.35), (-.6, 0.28),
                             connectionstyle="arc3,rad=-.2", **kw)
c1 = patches.FancyArrowPatch((.72, -0.27), (-.64, -0.43),
                             connectionstyle="arc3,rad=-.55", **kw)
#(tail coordinates),(head coordinates)

ax.text(0, 0, 'Tomato Curly\n Stunt Virus\n 2,766 bp', ha='center',
     bbox=dict(boxstyle='round', edgecolor='none', color='white'),
     )

plt.setp(bases+ inside + outside + bp_ticks, width=width)

#add patches of arrows to graph
for x in [v1,v2,c3,c2,c1]:
    plt.gca().add_patch(x)

#legend
gene_names = ['Coat Protein V1', 'V2 Protein', 'C3 Protein','C2 Protein','C1 Protein']

custom_lines = [Line2D([0], [0], color='slateblue', lw=4),
                Line2D([0], [0], color='darkblue', lw=4),
                Line2D([0], [0], color='cornflowerblue', lw=4),
                Line2D([0], [0], color='blueviolet', lw=4),
                Line2D([0], [0], color='blue', lw=4)]

plt.legend(custom_lines, gene_names, loc="upper right", bbox_to_anchor=(1.3,1.1), prop={'size': 7})



#save as png file
plt.savefig("tomato_curly_stunt_virus_genome.png") 

plt.show()



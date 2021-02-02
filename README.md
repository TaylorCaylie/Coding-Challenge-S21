# ACM Research Coding Challenge (Spring 2021)

## No Collaboration Policy

**You may not collaborate with anyone on this challenge.** You _are_ allowed to use Internet documentation. If you _do_ use existing code (either from Github, Stack Overflow, or other sources), **please cite your sources in the README**.

## Submission Procedure

Please follow the below instructions on how to submit your answers.

1. Create a **public** fork of this repo and name it `ACM-Research-Coding-Challenge-S21`. To fork this repo, click the button on the top right and click the "Fork" button.
2. Clone the fork of the repo to your computer using `git clone [the URL of your clone]`. You may need to install Git for this (Google it).
3. Complete the Challenge based on the instructions below.
4. Submit your solution by filling out this [form](https://acmutd.typeform.com/to/uqAJNXUe).

## Question One

Genome analysis is the identification of genomic features such as gene expression or DNA sequences in an individual's genetic makeup. A genbank file (.gb) format contains information about an individual's DNA sequence. The following dataset in `Genome.gb` contains a complete genome sequence of Tomato Curly Stunt Virus. 

**With this file, create a circular genome map and output it as a JPG/PNG/JPEG format.** We're not looking for any complex maps, just be sure to highlight the features and their labels.

**You may use any programming language you feel most comfortable. We recommend Python because it is the easiest to implement. You're allowed to use any library you want to implement this**, just document which ones you used in this README file. Try to complete this as soon as possible.

Regardless if you can or cannot answer the question, provide a short explanation of how you got your solution or how you think it can be solved in your README.md file. However, we highly recommend giving the challenge a try, you just might learn something new!

## Circular Genome Representation 

There are two types of genome mapping - physical and genetic. I chose physical mapping based on the contents of the Genbank file and utilized the base span of the CDS to visualize the regions of nucleotides that correspond to the sequence of amino acids in a protein product measured by base pairs. 

I first opened the Genbank file for parsing using the Biopython package and the Sequence Input/Output interface. From the file the base span and sequence is obtained and stored. The virus being analyzed has a length of 2766 bp (base pairs) and 5 known gene sequences.

To plot this information in a circular format I used matplotlib to map out the protein products and their associated base spans. I did a multilevel graph to show overlapping protein products and used patches to visualize if the protein was on a complimentary strand. The legend displays a color coded guide to show the correlated proteins. 

The most outer circular layers shows the genome sequence. Adenine is blue, Cytosine is red, Guanine is green, and Thymine is yellow. The second most outer layer is a ruler using bp as a measurement to evaluate the corresponding sequences and protein products.

![alt text](https://github.com/TaylorCaylie/Coding-Challenge-S21/blob/main/tomato_curly_stunt_virus_genome.png)

## Resources Used

- Biopython was used to access the genbank file and obtain the cds and associated attributes. 
http://biopython.org/DIST/docs/biosql/python_biosql_basic.html
- Matplotlib allowed for the patches to visualize the direction and the plot itself.  
J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, 2007.

===============================================================================
         TarHunter: a tool for predicting conserved miRNA target 
                    and target mimics in plants.
===============================================================================

Version 1.0

Author: Xuan Ma <skyxma@tjnu.edu.cn>
 
This program is free software: you can redistribute it and/or modify it under the 
terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. See the GNU General Public License for more details.




===============================================================================
                                  Description
===============================================================================
TarHunter's ortho_mode is designed for the prediction of conserved CDS targets. Its 
strategy is that orthologous miRNAs target the corresponding positions of multiple 
sequence alignment of orthologous genes. It is performed in the following four steps:

(i)   orthologous miRNA clustering
(ii)  orthologous gene clustering
(iii) miRNA target searching by FASTA/RNAhybrid
(iv)  cross-species conservation analysis of target sites

TarHunter's homo_mode is designed for the prediction of conserved CDS/eTM targets.
It requires that 50-bp upstream and downstream of the target sites are cross-species 
conserved. It is performed in three steps: (i), (iii) and (iv) as mentioned above. 
eTM prediction is based on the criteria that the central mispairs are flanked by two 
highly paired regions.

Citation: Ma X. et al. TarHunter, a tool for predicting conserved microRNA targets 
and target mimics in plants. Bioinformatics. 2017 Dec. 11 


===============================================================================
                                 Prerequisites
===============================================================================
OS Platform: Linux
Perl 5.10, or above

Executables under PATH:
CDHIT        (http://weizhongli-lab.org/cd-hit/)
FASTA        (fasta36, http://faculty.virginia.edu/wrpearson/fasta/fasta36/)
GNU Parallel (http://www.gnu.org/software/parallel/)
MCL          (http://micans.org/mcl/)
RNAhybrid    (http://bibiserv.techfak.uni-bielefeld.de/rnahybrid)




===============================================================================
                                  Simple test
===============================================================================
Run the following commands:
cd test/

    perl ../TarHunter.pl -q mir.txt -s spe1.txt -b cds_db/ -o out




===============================================================================
                                    Examples
===============================================================================
cd test/

A. CDS target prediction (ortho_mode, one step)

    perl ../TarHunter.pl -q mir.txt -s spe1.txt -b cds_db/ -o out1


B. CDS target prediction (ortho_mode, two steps, recommended)

   B1. Produce orthologous gene alignment file
   
       perl ../TarHunter.pl -q mir.txt -s spe1.txt -b cds_db/ -p 0 -o ortho
       
   B2. Search conserved miRNA targets
   
       perl ../TarHunter.pl -q mir.txt -s spe1.txt -b cds_db/ -a ortho/ortho_aln_MUSCLE.afa -o out2


C. CDS target prediction (ortho_mode) using varied parameters

   (Note: it is recommended to run B1 before running the following commands, and add -a ortho/ortho_aln_MUSCLE.afa parameter)

   C1. Set score cutoff to 5
   
       perl ../TarHunter.pl -q mir.txt -s spe1.txt -b cds_db/ -f 5 -o out3
      (note: C1 captures miR172-MYB84 that is not detected in B2, as the default score cutoff is 4)
  
   C2. Use total mispair cutoff
   
       perl ../TarHunter.pl -q mir.txt -s spe1.txt -b cds_db/ -M 4 -o out4

   C3. Use RNAhybrid
   
       perl ../TarHunter.pl -q mir.txt -s spe1.txt -b cds_db/ -p 2 -o out5
    
   C4. Without conservation filter
   
       perl ../TarHunter.pl -q mir.txt -s spe1.txt -b cds_db/ -N -o out6
    
   C5. Without orthologous miRNA search
   
       perl ../TarHunter.pl -q mir.txt -s spe1.txt -b cds_db/ -R -o out7

       
D. CDS target prediction (homo_mode) using varied parameters
   
   D1. Search conserved miRNA targets
   
       perl ../TarHunter.pl -q mir.txt  -s spe1.txt  -b cds_db/ -G -o out8
       
   D2. Set identity cutoff to 50%
       
       perl ../TarHunter.pl -q mir.txt  -s spe1.txt  -b cds_db/ -G -c 0.5 -o out9
       (note: a loose identity cutoff (-c 0.5) groups AP2,TOE,SMZ together)
   
E. target mimics prediction

    perl ../TarHunter.pl -q mir.txt -s spe2.txt -b nc_db/ -I -o out10

   
F. noncoding target prediction (homo_mode)

    perl ../TarHunter.pl -q mir.txt  -s spe3.txt  -b nc_db/ -G -f 5 -c 0.5 -o out11
    (note: F use moderately relaxed settings (-f 5 -c 0.5) to capture miR390-TAS3)

   
G. prediction of the targets of non-miRBase miRNAs

    perl ../TarHunter.pl -q mir.txt -n ../miRs/sit_mir.fa -s spe4.txt -b cds_db/ -o out12
    (note: G predicts the targets of S. italica miRNAs (Chavez et al 2014 Nat Commun) that are not collected in miRBase)



   
===============================================================================
                                Important notes
===============================================================================
1. The names of the sequence files should be species name abbreviation followed
   by dot fa, e.g., ath.fa.
   
2. miRNA list can either contain miRNA IDs from multiple species, or contain 
   miRNA IDs from a single species (TarHunter can automatically search for 
   orthologous miRNAs in other species.
   
3. species name abbreviations and miRNA IDs should follow miRBase nomenclature.
   
4. Make sure the TarHunter script, the miRs and bin folders are under the same directory.

5. Coding sequences should be in +1 frame, otherwise will be neglected.
   
6. It is recommended to use FASTA (default search engine) in target searching, 
   as FASTA is much faster than RNAhybrid.

7. The script has been tested on a Fedora (release 16) server with Perl v5.14.3,
   and on a UBUNTU (release 14.04) PC with Perl v5.18.2.


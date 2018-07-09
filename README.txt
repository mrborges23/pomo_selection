# pomo_selection
This algorithms was built to quantify the patters of GC-bias gene conversion is great apes’ populations. Our strategy consisted in estimating the population parameter (mutation rates and selection coefficietns) using the stationary distribution of a multivariate Moran model with boundary mutations and allelic selection. The input data is the vector of allele counts, you can check how to buil in point (3). The inferential framework is use is Bayesian. 

1. Before running the algortihms or the examples provided please make sure you have the R package “MCMCpack” installed!

2. You can follow the example.R code for some simulated and real examples.

3. How to build the allele counts vector.
Consier an alignment of 20 sites for 4 individuals of the same population.

      1       5   * * * 10    * * 15    *   20
Ind1  C G T A A G G G T C G T C T G C T A T A
Ind2  C G T A A G A G G C G T C T G C T C T A
Ind3  C G T A A G A C T C G T T T G C T C T A
Ind4  C G T A A G A C T C G T C A G C T C T A

We list all the possible states, considering the number of individuals in study, and build the count vector by simply counting the number of times a certain state appears in the alignment. 

States	Counts
{A}        	3 (sites 4, 5 and 20) 		
{C}        	3 (sites 1, 10 and 16)
{G}		4 (sites 2, 6, 11 and 15) 
{T}		4 (sites 3, 12, 17 and 19)				
{1A:3C}	1 (site 18)
{2A:2C}	0
{3A:1C}	0
{1A:3G}	0
{2A:2G}	0
{3A:1G}	1 (site 7)
{1A:3T}	1 (site 14)
{2A:2T}	0
{3A:1T}	0
{1C:3G}	0
{2C:2G}	1 (site 8)
{3C:1G}	0
{1C:3T}	0
{2C:2T}	0
{3C:1T}	1 (site 13)
{1G:3T}	1 (site 9)
{2G:2T}	0
{3G:1T}	0

Note: The order of the monomorphic {A,C,G,T} and the polymorphic {1A:(N-1)C,...(N-1)A:1C,1A:(N-1)G,...,(N-1)A:1G,1A:(N-1)T,…,(N-1)A:1T,1C:(N-1)G,...,(N-1C):1G,1C:(N-1)T,...,(N-1)C:1T,1G:(N-1)T,…,(N-1)G:1T} states should be respected; incorrect labels may be returned if this order is not respected!

Check the input file counts_pongo_pygmaeus.txt. 

4. Cite us:


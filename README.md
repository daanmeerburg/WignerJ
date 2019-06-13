# WignerJ
fast algebraic wignerJ for (l1 l2 l3 0 0 0)

Just run ifrt or gfortran. No fancy stuff. 

Compares widely used itterative verion to simple algebraic version which is based on a better 
approximation of the Sterling approximation 
We find > 1% agreememt (more like .5%) on all triplets. 

There are some notes here:

https://www.overleaf.com/read/hsymsttscxxr

The code will run and compare the two methods and prints triplets, results and relative error in %. 
You can increase imax if you want to check higher values (they tend to get better)


June 13, 2019

Added algebraic wignerJ for (l1 l2 l3 -2 0 2) based on Eq. B1 from arXiv 0001303
Results are not that good for squeezed triplets and near zero crosses. However, majority of (odd!) triplets
leads to <3% errors. Nice thing is that effectively it is just multiplicative of the other WignerJ so very
little extra code.  

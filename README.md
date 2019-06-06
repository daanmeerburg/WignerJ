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


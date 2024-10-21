Sharkovski ordering in a delayed perturbation of the Roessler ODE
=================================================================

This a program to accompany paper 

> Sharkovskii theorem for Delay DifferentialEquations and other 
> infinite dimensionaldynamical systems.
>
> by Anna Gierzkiewicz and Robert Szczelina (2024)
> [Preprint]

If you have the ```doi``` version of this repository, 
check the github 
[https://github.com/robsontpm/sharkovski-roessler-dde]
for eventual changes.


Compilation
-----------

The project requires downloading and installing capdDDEs library: 
[https://github.com/robsontpm/capdDDEs]. The library itself uses
CAPD [http://capd.ii.uj.edu.pl/], but it is included in the repo
already. 

Steps to compile and run programs:

```
# clone the parent repository
git close https://github.com/robsontpm/capdDDEs

# build that repository
cd capdDDEs/
bash tldr.sh

# go to the right place to store your own programs
# this will make use of the capdDDEs compilation infrastructure
cd programs/results

# clone this repository
git clone https://github.com/robsontpm/sharkovski-roessler-dde
cd sharkovski-roessler-dde

# build all programs
make all

# run the proof
cd scripts
bash RUN_PROOF.sh
```

NOTE: The running time is approximately 3h on ~3GHz PC with 10 availabe 
physical cores. If you have more cores, you can edit ```RUN_PROOF.sh```
to change ```NUM_JOBS=10``` into ```NUM_JOBS=20``` (if you have 20 cores).


Troubleshooting
---------------

The codes will most likely not compile on ```Windows```. 

The codes will most likely not work on ```Mac```. hey should
compile, but you might need to generate data from scratch. 
See below for detials.  

The proof was prepared and run on ```Ubuntu 20.04.6 LTS```, 
compiled with ```gcc (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0```.
Any PC with x86 architecture should be able to run the
proof. In a rare case, that the proof fails, it might be
needed to recreate data from scratch. 
See ```README.md``` in ```scripts/``` to see how to do that.

In case of any problems, do not hesitate to contact 
us by email: robert(dot)szczelina(at)uj(dot)edu(dot)pl. 


Project structure
-----------------

- ```bin```
  this is where compiled programs go
- ```scripts```
  this is where the protocol of the proof is established
  as a collection of bash/python scripts.
  Since proof is more involved than the version for ODEs,
  especially the generating the data for proofs, we
  decided to automate it for user convenience. We use
  scripting languages as they are easier to read and follow,
  brake down into parts that can be run separatley
  than a very complicated one-file-does-everything program in C++.
  For the sake of completeness and reproductibility, we include
  the full protocol described in the manuscript, so that 
  anyone can recreate the data from scratch or use/modify it
  to do the proof for other interesting cases. 
- ```proof-data```
  data and results from the following the protocol on our 
  computers. You can check the proof on yours computer just
  by doing the last step - veryfication
- ```proof-output```
  here the output of the proof will be stored. For convenience,
  we put here the data obtained on our PC. Also, there are some 
  ```gnuplot``` scripts to generate nice pictures used in the 
  accompanying paper.
- ```tmp``` this directory will show up only if you are
  regenerating data from scratch, following the protocol
  in ```scripts```. 
- ```Makefile```
  This guides the compilation. Its structure is such
  that it works if one put the repository 
  inside the programs/results of the main ddesCAPD package. 
- ```.h``` and ```.cpp``` files - the main ingredients of 
  the computer assisted part. In header files are common 
  subroutines, and source code files should have names
  describing their use. Each source file, except of ```utils.cpp```
  will produce separate program in the ```/bin``` directory. 
  See how they are used in ```proof-scripts```
- ```.gitignore``` this is a standard ```git``` configuration 
  file
- ```README.md```
  this file, in markup format, that can be pretty-formated by
  various programs (github, sublime text, etc).
- ```TODO.md```
  developers comments
  



Sharkovski ordering in a delayed perturbation of the Roessler ODE
=================================================================

***FOR DEVS!***: Please clone this into the ```programs/results```
or ```programs/devel``` directory in the capdDDEs main repository!


Project structure
-----------------

- ```bin```
  this is where compiled programs go
- ```proof-scripts```
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
  



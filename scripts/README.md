The protocol
============


1. The ```RUN_GENCOORS.sh``` is the first to be run:

   ```
   bash RUN_GENCOORDS.sh
   ```

   It will create ```../tmp/``` directory. Inside, a bunch of
   ```.vector.bin``` and ```imatrix.bin``` files will be. Those
   are the coordinate frame for the main set - the
   trapping region for the equation. It will also generate some
   ```.dat``` files to be used for pictures. Those will be automatically 
   collected and used in subsequent phases. 

2. ```bash RUN_INV_M.sh```
   This will compute rigorously the matrix $M^{-1}$.
   The process takes several minutes. 

3. ```bash RUN_PREP.sh```
   This will try to find initial set. The run is long
   (about 3h on 10 cores). 
   Run ```bash RUN_PREP.sh``` second time to check
   and to collect data. Again, 3h. 
   In fact, if this step is successful, you have proved
   the existence of isolating neighbourhood. 

4. ```bash RUN_PROOF.sh```
   this will run the final proof. Again, ~3h. 


Other files and configuration
-----------------------------

The other scripts are as follows:
- ```python3 GEN_JOBS.py``` (needs python3)
  if you change the configuration inside this file
  you need to run this. Or if you want to increase or 
  decrease the number of pieces into which divide the 
  sets. After that you essentially need to rerun everything
  from the list above. 

- ```PREP_JOBLIST.sh``` and ```PROOF_JOBLIST.sh```
  those files are created by ```python3 GEN_JOBS.py```. 
  You can copy single lines from those files and run 
  them separately, e.g. when debugging your own proofs. 

- ```GEN_RUN_HULL.sh```
  if you change the number of divisions of the 
  isolating neighbourhood G, then you need to run this file
  once in step 3 of the protocol, after first invocation
  of ```bash RUN_PREP.sh```. Then you need to redo step 3,
  this time you do not need to invoke ```GEN_RUN_HULL.sh```
  You can experiment with the context of this file,
  especially the scaling factors. 

- ```RUN_HULL.sh```
  generated with ```GEN_RUN_HULL.sh```, it is used 
  in step 3 in ```bash RUN_PREP.sh``` to produce
  initial estimates on the tail (far and close) of the C^n_p 
  representations.



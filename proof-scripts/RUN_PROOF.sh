NUM_JOBS=10
cd bin
cat ../JOBLIST.txt | xargs -P $NUM_JOBS -I{} bash -c "{} 2> /dev/null"

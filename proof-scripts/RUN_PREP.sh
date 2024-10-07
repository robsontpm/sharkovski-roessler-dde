NUM_JOBS=10
cd bin
cat ../PREPLIST.txt | xargs -P $NUM_JOBS -I{} bash -c "{} 2> /dev/null"

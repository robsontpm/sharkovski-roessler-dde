# you can change this to the number of computational cores you have
NUM_JOBS=10
# this is where output data from proof will be (should be the same as those from PROOF_JOBLIST.sh)
OUTDIR="../proof-output"
# this is where the log will be
OUTFILE="$OUTDIR/computation.log"

# make sure the file is empty, put date of the run in it.
date > $OUTFILE
# we run jobs in parallel, using xargs
# this is not so great, bbut for now we have not figured out
# how to use OpenMP for this. I have tried, but failed
# with the capdDDEs code. The basic capd seem to work ok.
cat ./PROOF_JOBLIST.sh | xargs -P $NUM_JOBS -I{} bash -c "{} 2> /dev/null" | tee -a $OUTFILE

# check if all the conditions are satisfied
if grep -q "false" "$OUTFILE"; then
  echo "FAILURE: one of the conditions not meet."  
else
  echo "SUCCESS: all conditions of the theorem confirmed."  
fi
# we always show the user where the log is
echo "Check log in '$OUTFILE'"


# gather data for pictures:
find $OUTDIR/ -mindepth 1 -type d | xargs -I{} bash -c "echo "" > {}-all.dat"
find $OUTDIR/ -mindepth 1 -type d | xargs -I{} bash -c "find {} -type f | grep '.dat' | xargs -I XX grep -v '^\s*#' XX >> {}-all.dat"


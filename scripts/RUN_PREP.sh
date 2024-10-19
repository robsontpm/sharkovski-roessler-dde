NUM_JOBS=10

WD="../tmp"
# this is where output data from proof will be (should be the same as those from PROOF_JOBLIST.sh)
OUTDIR="../proof-data"

if [ -f $WD/KEEP_LOOKING.txt ]; then
    rm $WD/KEEP_LOOKING.txt
fi

cat ./PREP_JOBLIST.sh | xargs -P $NUM_JOBS -I{} bash -c "{} 2> /dev/null"

if [ -f $WD/KEEP_LOOKING.txt ]; then
    bash RUN_HULL.sh
    echo "KEEP LOOKING"
else
    echo "DONE!"
    cp $WD/C_0_x0.ivector.bin $OUTDIR
    cp $WD/C_1_x0.ivector.bin $OUTDIR
    cp $WD/C_2_x0.ivector.bin $OUTDIR
    cp $WD/G_invM.imatrix.bin $OUTDIR
    cp $WD/G_M.imatrix.bin $OUTDIR
    cp $WD/G_r0.ivector.bin $OUTDIR
    cp $WD/G_x0.ivector.bin $OUTDIR
    cp $WD/G_Xi.ivector.bin $OUTDIR
fi
# where the helper scrip convmatrix is
SCRIPTDIR=../../../../bin
# where to find the sets
SRCDIR=../tmp/estimate_piece
# where to put result
DSTDIR=../tmp
# should be setup as in the programs! M and d*p, respectively!
DIM_R0=771
DIM_XI=192
# current date for backup


# code to make backup, and also clear the file
echo 'DATE=`date +"%Y-%m-%d--%H-%M"`' > RUN_HULL.sh 
echo "cp $DSTDIR/G_r0.ivector.bin $DSTDIR/bak-\$DATE--G_r0.ivector.bin" >> RUN_HULL.sh 
echo "cp $DSTDIR/G_r0.ivector.bin $DSTDIR/bak-\$DATE--G_Xi.ivector.bin" >> RUN_HULL.sh 
echo "" >> RUN_HULL.sh

# this computes hull of r0 vectors
echo "$SCRIPTDIR/vectorhull ${DIM_R0} bin $DSTDIR/G_r0.ivector.bin \\" >> RUN_HULL.sh
find $SRCDIR -type f | grep new_r0 | xargs -I{} echo "    {} \\" >> RUN_HULL.sh
echo "" >> RUN_HULL.sh

# this computes hull of Xi vectors
echo "$SCRIPTDIR/vectorhull ${DIM_XI} bin $DSTDIR/G_Xi.ivector.bin \\" >> RUN_HULL.sh
find $SRCDIR -type f | grep new_Xi0 | xargs -I{} echo "    {} \\" >> RUN_HULL.sh
echo "" >> RUN_HULL.sh

# this grows the vectors a little bit
echo $SCRIPTDIR/mulvector ${DIM_R0} bin $DSTDIR/G_r0.ivector.bin '[-1.5,1.5]' $DSTDIR/G_r0.ivector.bin >> RUN_HULL.sh
echo $SCRIPTDIR/growvector ${DIM_XI} bin $DSTDIR/G_Xi.ivector.bin '1.5' $DSTDIR/G_Xi.ivector.bin >> RUN_HULL.sh 

# this produces human readable vectors
echo $SCRIPTDIR/convvector ${DIM_R0} bin $DSTDIR/G_r0.ivector.bin conv interval $DSTDIR/G_r0.ivector >> RUN_HULL.sh
echo $SCRIPTDIR/convvector ${DIM_XI} bin $DSTDIR/G_Xi.ivector.bin conv interval $DSTDIR/G_Xi.ivector >> RUN_HULL.sh 
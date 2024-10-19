# where find gencoords.log
WD="../tmp"
# where to put output
OUTDIR="../proof-output"
# this is where the log will be
LOGFILE="$WD/gencoords.log"
# this is where the log will be
OUTFILE="$OUTDIR/C_coords.gp"

../bin/Roessler_DDE_gencoords "wd=$OUTDIR/" verbose 2>&1 | tee $LOGFILE

echo "" > $OUTFILE

echo '# those are the coordiates translations needed to move the resulting box' >> $OUTFILE
echo '# to the respective center. This must be here, because of the way the program' >> $OUTFILE
echo '# now returns the images in good coordinates (around those points)' >> $OUTFILE
echo '# So, the sets are 0,0- centered on output.' >> $OUTFILE
echo '# To put them in good perspective of the full set, we need to translate them.' >> $OUTFILE

for i in 0 1 2
do
    grep "mid in coordinates" $LOGFILE \
        | grep "C\[$i\]" \
        | grep -o "\{\[[^\ ]*," \
        | grep -o "[^\{\[,]*" \
        | xargs -I{} echo "C_${i}_X = {}" >> $OUTFILE
    grep "mid in coordinates" $LOGFILE \
        | grep "C\[$i\]" \
        | grep -o ",\[[^\ ]*," \
        | grep -o "[^\{\[,]*" \
        | xargs -I{} echo "C_${i}_Y = {}" >> $OUTFILE
done

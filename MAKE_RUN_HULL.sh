echo "" > RUN_HULL.sh 

DIM_R0=771
DIM_XI=192
echo "../../../bin/vectorhull ${DIM_R0} bin bin/grid3_new_r0.ivector.bin \\" >> RUN_HULL.sh
find ./bin/estimate_piece -type f | grep new_r0 | xargs -I{} echo "{} \\" >> RUN_HULL.sh
echo "" >> RUN_HULL.sh

echo "../../../bin/vectorhull ${DIM_XI} bin bin/grid3_new_Xi.ivector.bin \\" >> RUN_HULL.sh
find ./bin/estimate_piece -type f | grep new_Xi0 | xargs -I{} echo "{} \\" >> RUN_HULL.sh
echo "" >> RUN_HULL.sh

echo ../../../bin/convvector ${DIM_R0} bin bin/grid3_new_r0.ivector.bin conv interval bin/grid3_new_r0.ivector >> RUN_HULL.sh
echo ../../../bin/convvector ${DIM_XI} bin bin/grid3_new_Xi.ivector.bin conv interval bin/grid3_new_Xi.ivector >> RUN_HULL.sh 

echo ../../../bin/mulvector ${DIM_R0} bin bin/grid3_new_r0.ivector.bin '[-2,2]' bin/grid3_new_r0.ivector.bin >> RUN_HULL.sh
echo ../../../bin/growvector ${DIM_XI} bin bin/grid3_new_Xi.ivector.bin '2.0' bin/grid3_new_Xi.ivector.bin >> RUN_HULL.sh 

echo ../../../bin/convvector ${DIM_R0} bin bin/grid3_new_r0.ivector.bin conv interval bin/grid3_new_r0.ivector >> RUN_HULL.sh
echo ../../../bin/convvector ${DIM_XI} bin bin/grid3_new_Xi.ivector.bin conv interval bin/grid3_new_Xi.ivector >> RUN_HULL.sh 
find ./bin/Pimages/ -type f | grep ".dat" | grep Pc3_0 | xargs -I{} bash -c "cat {} && echo ''" > bin/plots/Pc3_0.dat
find ./bin/Pimages/ -type f | grep ".dat" | grep Pc3_1 | xargs -I{} bash -c "cat {} && echo ''" > bin/plots/Pc3_1.dat
find ./bin/Pimages/ -type f | grep ".dat" | grep Pc3_2 | xargs -I{} bash -c "cat {} && echo ''" > bin/plots/Pc3_2.dat
find ./bin/Pimages/ -type f | grep ".dat" | grep grid3 | xargs -I{} bash -c "cat {} && echo ''" > bin/plots/Pgrid3.dat
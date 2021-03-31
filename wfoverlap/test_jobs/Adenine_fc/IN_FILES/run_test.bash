$OVDIR/../bin/wfoverlap.x -f ciov_all.in > ciov_all.out
$OVDIR/../bin/wfoverlap.x -f ciov_fc.in  > ciov_fc.out
$OVDIR/../bin/wfoverlap.x -f dyson.in > dyson.out

export OMP_NUM_THREADS=2
$OVDIR/../bin/wfoverlap.x -f ciov_fc.in -m 7  > ciov_low-mem.out

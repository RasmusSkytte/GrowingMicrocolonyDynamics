# Compiler variables

all :
	make process
	make videos
	make analyze

process :

	### Process the raw experiment images ###
	#
	#matlab -nodesktop -nosplash -r "run('../data/ProcessExportedData.m'); exit" > /dev/null 2>&1
	#
	find experiments -name '*.tif' | xargs -i convert {} -compress lzw {}
	#

LIST=B1X1 B1X2 B1X3 B1X4 B1X5 B2X1 B2X2 B2X3 B2X4 B2X5 B3X1 B3X2 B3X3 B3X4 B3X5 B4X1 B4X2 B4X3 B4X4 B4X5 B5X1 B5X2 B5X3 B5X4 B5X5 B6X1 B6X2 B6X3 B6X4 B6X5
videos :
	for p in  $(LIST); \
	do \
		echo $$p | xargs -i ffmpeg -framerate 5  -i 'experiments/Processed data/{}_t%03d.tif' -s:v 1280x720 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p experiments/Videos/{}.mp4 -y ; \
	done

analyze :
	# matlab -nodesktop -nosplash -r "run('analysis/functions/minimizerGrowthCurve.m'); exit"  > Log_MinimizeGrowthCurve.txt &
	# matlab -nodesktop -nosplash -r "run('analysis/functions/determineGrowthCurveUncertainties.m'); exit"
	matlab -nodesktop -nosplash -r "run('analysis/functions/minimizerPhageAttack.m'); exit" > Log_MinimizePhageAttack.txt &
	matlab -nodesktop -nosplash -r "run('analysis/functions/minimizerPhageAttackCoupled.m'); exit" > Log_MinimizePhageAttackCoupled.txt &
	matlab -nodesktop -nosplash -r "run('analysis/functions/minimizerPhageAttackSingle.m'); exit" > Log_MinimizePhageAttackSingle.txt &
	#exit" 2>&1 | tee DataSet4_log.txt

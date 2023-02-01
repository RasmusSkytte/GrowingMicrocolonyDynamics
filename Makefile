all :
	make process
	make videos
	make filter
	make growth_curve
	make transparency_parameter
	make phage_attck
	make visualize

process :

	### Process the raw experiment images ###
	matlab -nodesktop -nosplash -r "run('../data/ProcessExportedData.m'); exit" > /dev/null 2>&1
	find experiments -name '*.tif' | xargs -i convert {} -compress lzw {}

LIST=B1X1 B1X2 B1X3 B1X4 B1X5 B2X1 B2X2 B2X3 B2X4 B2X5 B3X1 B3X2 B3X3 B3X4 B3X5 B4X1 B4X2 B4X3 B4X4 B4X5 B5X1 B5X2 B5X3 B5X4 B5X5 B6X1 B6X2 B6X3 B6X4 B6X5
videos :
	for p in  $(LIST); \
	do \
		echo $$p;\
		echo $$p | xargs -i ffmpeg -framerate 5  -i 'experiments/Processed data/{}_t%03d.tif' -s:v 1280x720 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p experiments/Videos/{}.mp4 -y ; \
	done

filter :
	matlab -nodesktop -nosplash -r "run('analysis/functions/filterData.m'); exit"

growth_curve :
	matlab -nodesktop -nosplash -r "run('analysis/functions/minimizerGrowthCurve.m'); exit"
	matlab -nodesktop -nosplash -r "run('analysis/functions/determineGrowthCurveUncertainties.m'); exit"

transparency_parameter :
	matlab -nodesktop -nosplash -r "run('analysis/functions/minimizerTransparencyParameterAndThresholds.m'); exit"

phage_attck :
	matlab -nodesktop -nosplash -r "run('analysis/functions/minimizerPhageAttack.m'); exit" 
	matlab -nodesktop -nosplash -r "run('analysis/functions/testPhageAttack.m'); exit" 

visualize :
	matlab -nodesktop -nosplash -r "run('analysis/functions/testGrowthCurve.m'); exit" 
	matlab -nodesktop -nosplash -r "run('analysis/functions/testPhageAttack.m'); exit" 
	matlab -nodesktop -nosplash -r "run('analysis/functions/testBalance.m'); exit" 
	matlab -nodesktop -nosplash -r "run('analysis/functions/testBalanceSensitivity.m'); exit" 
	matlab -nodesktop -nosplash -r "run('analysis/functions/testTransparencyParameterAndThresholds.m'); exit" 

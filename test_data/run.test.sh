source /nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/common/Anaconda3/bin/activate env_nf

export NXF_SINGULARITY_CACHEDIR=/nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/test/ASEprep/singularity_images

nextflow run main.nf \
        -profile singularity,gl \
        -params-file test_data/params.test.yaml \
	-resume \
        #-with-dag test_data/reports_test/run_report/dag.png \
	#-with-trace test_data/reports_test/run_report/trace.txt \
	#-with-report test_data/reports_test/run_report/report.html \
	#-with-timeline test_data/reports_test/run_report/timeline.html \

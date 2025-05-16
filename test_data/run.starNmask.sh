source /nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/common/Anaconda3/bin/activate env_nf

export NXF_SINGULARITY_CACHEDIR=/nfs/mm-isilon/bioinfcore/ActiveProjects/weishwu/test/ASEprep/singularity_images

nextflow run main.nf \
        -profile lh,singularity \
	-params-file test_data/params.starNmask.yaml \
	-resume \
        -with-dag test_data/reports_starNmask/run_report/dag.png \
	-with-trace test_data/reports_starNmask/run_report/trace.txt \
	-with-report test_data/reports_starNmask/run_report/report.html \
	-with-timeline test_data/reports_starNmask/run_report/timeline.html \
	#-with-dag test_data/reports/run_report/dag.png \
	#-resume \
	#-preview

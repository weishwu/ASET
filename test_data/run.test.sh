export NXF_SINGULARITY_CACHEDIR=~/sifs/

nextflow run main.nf \
	-c config/test.config \
        -profile singularity,cluster \
        -params-file test_data/params.test.yaml \
	-resume \
        -with-dag test_data/reports_test/run_report/dag.png \
        -with-trace test_data/reports_test/run_report/trace.txt \
        -with-report test_data/reports_test/run_report/report.html \
        -with-timeline test_data/reports_test/run_report/timeline.html \

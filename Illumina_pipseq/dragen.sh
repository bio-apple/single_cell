/opt/dragen/4.4.6/bin/dragen --fastq-list ${1} \
--ref-dir=${2} --annotation-file=${3} \
--output-directory ${4} \
--output-file-prefix ${5} \
--scrna-enable-pipseq-mode=true \
--fastq-list-sample-id ${5}

/usr/bin/dragen-reports \
-f -d ${4} -o ${4}/report.html \
-m /opt/dragen-reports/manifests/scrna.json
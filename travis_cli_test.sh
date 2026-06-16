MIGMAP="java -Xmx4G -jar `ls build/libs/migmap-*.jar` --blast-dir . --data-dir data/ -S human -R IGH --details all --all-alleles"
$MIGMAP post/raji_R12.fastq.gz raji.txt
$MIGMAP --by-read post/raji_R12.fastq.gz - > raji_reads.txt

cat raji.txt | head

MIGMAP_MODULE="java -Xmx4G -cp `ls build/libs/migmap-*.jar`"

$MIGMAP_MODULE com.antigenomics.migmap.MergeContigs raji.txt raji_mc.txt
$MIGMAP_MODULE com.antigenomics.migmap.Correct raji_mc.txt raji_mc_corr.txt
$MIGMAP_MODULE com.antigenomics.migmap.Analyze -S human -R IGH raji_mc_corr.txt raji

ls -lh

for f in raji.txt raji_reads.txt raji_mc.txt raji_mc_corr.txt raji.edge.txt raji.net.txt raji.shm.txt
do
    if [[ ! -s $f ]]
        then exit 1
    fi
done
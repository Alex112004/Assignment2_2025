set -e
sed '1d' UTF-8assigment2.csv | while read meta
do
  accn=$(echo $meta | cut -d "," -f 1)
  name=$(echo $meta | cut -d "," -f 32)
  echo $name
  fasterq-dump ${accn}
  cat ${accn}_1.fastq | gzip > ${name}_1.fastq.gz
  rm ${accn}_1.fastq
  cat ${accn}_2.fastq | gzip > ${name}_2.fastq.gz
  rm ${accn}_2.fastq
done


module load python/anaconda3.6
module load sratoolkit
module load java
module load rsem/1.3.0
module load STAR/2.6.0
module load parallel

######## install aspera if not installed #########
## CSV files
CSV_file=$1
## number of threads/processors
NUM_OF_THREAD=$2
## cell type
CELL=$3


## activate conda environment
source activate /projects/e31039/softwares/pythonenvs/finalproject

## function for parallel fastq download
# input: SRRXXXXXXX num_of_thread
function LoadDump {
  echo "* prefetch $1"
  prefetch -a "$HOME/.aspera/bin/ascp|$HOME/.aspera/etc/asperaweb_id_dsa.openssh" -O ./ -X 999999999 $1

  if [[ -e ${1}/${1}.sra ]]; then
    echo "* running fastq-dump in parallel ..."
    parallel-fastq-dump -s ${1}/${1}.sra -t $2 -O ./ --tmpdir ./ --split-3 && rm -R ${1}
  else
    echo '* [ERROR]' $1 'apparently not successfully loaded' && exit 1
  fi
}; export -f LoadDump


# download
if ! type "efetch" > /dev/null; then
  print "* Please install E-utilities."
fi

# if input GSM list
# GSM_LIST=$1
# GSMs=`cat $GSM_LIST|cut -f1`

# if input whole csv file
GSMs=$(cat $CSV_file | cut -d ',' -f 4 | sed 1d)

# create output directory
mkdir -p ./$CELL; cd ./$CELL

for GSM in $GSMs; do
  echo "* $GSM retrieved from NCBI GEO ..."
  if [ ! -d $GSM ]; then
    mkdir -p $GSM; cd $GSM
  else
    echo "* directory for $GSM already exists!"
    cd $GSM
    # if [[ $(find . -maxdepth 1 -type f -name "*.results" 2>/dev/null | wc -l) != 0]]; then
    if compgen -G "./*.results" > /dev/null; then
      echo "* RSEM already completed for: $GSM!"
      cd ..
      continue
    fi
  fi
  SRRs=`esearch -db sra -query $GSM |efetch -format docsum |xtract -pattern DocumentSummary -element Run@acc`
  num_SRRs=$(wc -w <<< $SRRs)
  if [[ $num_SRRs -eq "0" ]]; then
    echo "* ERROR: $GSM contains no SRRs! Skipped."
    cd ..
    continue # when no SRR included in GSM
  fi
  for SRR in $SRRs; do
    parallel --lb LoadDump ::: $SRR ::: $NUM_OF_THREAD
    echo "* fastq-dump completed!"
    echo " "
  done

  curr_dir=$(pwd)
  # determine single end or paired end, using compgen -G "<glob-pattern>" (https://stackoverflow.com/questions/2937407/test-whether-a-glob-has-any-matches-in-bash)
  if compgen -G "./*_1.fastq" > /dev/null; then # paired end
    if [[ $num_SRRs -gt "1" ]]; then # more than one SRR runs, merge
      echo "* Merging fastq files from individual runs..."
      for f in *_1.fastq ; do cat "$f" >> ${GSM}_1.fastq && rm "$f" || break ; done
      for f in *_2.fastq ; do cat "$f" >> ${GSM}_2.fastq && rm "$f" || break ; done
    else # rename
      mv *_1.fastq ${GSM}_1.fastq
      mv *_2.fastq ${GSM}_2.fastq
    fi
    echo "* Running trimmomatic..."
    java -jar /projects/e31039/softwares/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads $2 -trimlog ${curr_dir}/${GSM}_trimlog.txt ${curr_dir}/${GSM}_1.fastq ${curr_dir}/${GSM}_2.fastq -baseout ${curr_dir}/${GSM}.fastq ILLUMINACLIP:/projects/e31039/softwares/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    echo " "

    echo "* Running RSEM..."
    rsem-calculate-expression  -p $2  --no-bam-output --star --star-path /software/STAR/2.6.0/bin/Linux_x86_64_static --paired-end ${curr_dir}/${GSM}_1P.fastq ${curr_dir}/${GSM}_2P.fastq /projects/e31039/data/genomes_index/projectIdx/hg38StRsm ${curr_dir}/$GSM
    echo "* RSEM output saved to: $(pwd)/$GSM"
    echo " "

    echo "* Removing fastq and temp files..."
    rm -r ./${GSM%.bam}.temp ./${GSM%.bam}.stat
    rm *.fastq

    echo "* JOB FINISHED: $GSM"
    echo " "
  else # single end
    if [[ $num_SRRs -gt "1" ]]; then # more than one SRR runs, merge
      echo "* Merging fastq files from individual runs..."
      for f in *.fastq ; do cat "$f" >> ${GSM}.fastq && rm "$f" || break ; done
    else # rename
      mv *.fastq ${GSM}.fastq
    fi
    echo "* Running trimmomatic..."
    java -jar /projects/e31039/softwares/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads $2 -trimlog ${curr_dir}/${GSM}_trimlog.txt ${curr_dir}/${GSM}.fastq ${curr_dir}/${GSM}_P.fastq ILLUMINACLIP:/projects/e31039/softwares/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    echo " "

    echo "* Running RSEM..."
    rsem-calculate-expression --no-bam-output --star --star-path /software/STAR/2.6.0/bin/Linux_x86_64_static -p $2 ${curr_dir}/${GSM}_P.fastq /projects/e31039/data/genomes_index/projectIdx/hg38StRsm ${curr_dir}/$GSM
    echo "* RSEM output saved to: $(pwd)/$GSM"
    echo " "

    echo "* Removing fastq and temp files..."
    rm -r ./${GSM%.bam}.temp ./${GSM%.bam}.stat
    rm *.fastq
    rm *_trimlog.txt

    echo "* JOB FINISHED: $GSM"
    echo " "

  fi
  cd ..
done

#!/bin/bash    

for i in "$@"
do
case $i in
    -linc=*|--extension=*)
    linc_fasta="${i#*=}"
    shift # past argument=value
    ;;
    -pct=*|--searchpath=*)
    pcts_fasta="${i#*=}"
    shift # past argument=value
    ;;
    -filter=*|--lib=*)
    filter_answer="${i#*=}"
    shift # past argument=value
    ;;
    -pctGroup=*|--lib=*)
    pctGroup="${i#*=}"
    shift # past argument=value
    ;;
    --default)
    DEFAULT=YES
    shift # past argument with no value
    ;;
    *)
          # unknown option
    ;;
esac
done

lc(){
    case "$1" in
        [A-Z])
        n=$(printf "%d" "'$1")
        n=$((n+32))
        printf \\$(printf "%o" "$n")
        ;;
        *)
        printf "%s" "$1"
        ;;
    esac
}

f2i() {
  awk 'BEGIN{for (i=1; i<ARGC;i++)
   printf "%.0f\n", ARGV[i]}' "$@"
}

# Filter by biotype and transcript type
echo "Filter option: "
filter_answer= lc "$filter_answer"
echo ""

mkdir "temp"
mkdir "result"


 if [ "$filter_answer" == "y" ] || [ "$filter_answer" == "yes" ];
 then
   perl extractPctAndLinc.pl "$linc_fasta" "$pcts_fasta" "temp/lincs_filtered.fa" "temp/pcts_filtered.fa"
   linc_fasta="temp/lincs_filtered.fa"
   pcts_fasta="temp/pcts_filtered.fa"
 fi

# blast
mkdir blast
makeblastdb -in "$linc_fasta" -out "blast/lincs" -dbtype "nucl"
makeblastdb -in "$pcts_fasta" -out "blast/pcts" -dbtype "nucl"

blastn -query "$linc_fasta" -out "temp/blast_lincs" -outfmt 7 -evalue 1e-12  -max_target_seqs 2 -db "blast/pcts"
blastn -query "$pcts_fasta" -out "temp/blast_pcts" -outfmt 7 -evalue 1e-12  -max_target_seqs 2 -db "blast/lincs"

perl filterBlast.pl "temp/blast_lincs" "$linc_fasta" "temp/blast_filtered_lincs.fa"
perl filterBlast.pl "temp/blast_pcts" "$pcts_fasta" "temp/blast_filtered_pcts.fa"

linc_fasta="temp/blast_filtered_lincs.fa"
pcts_fasta="temp/blast_filtered_pcts.fa"

echo "BLAST FILTER processed"

# ORFS
./txCdsPredict -anyStart "$linc_fasta" "temp/lincs.cds"
./txCdsPredict -anyStart "$pcts_fasta" "temp/pcts.cds"
echo "ORFs features processed"
echo ""

# feature file
perl extractFeatures.pl "$linc_fasta" "temp/lincs.cds" "temp/linc_features.csv" "linc"
perl extractFeatures.pl "$pcts_fasta" "temp/pcts.cds"  "temp/pcts_features.csv" "pct"


# size of training and test gropu
size_lincs=$(wc -l "temp/linc_features.csv" | sed 's/[^0-9]*//g' )
size_pcts=$(wc -l "temp/pcts_features.csv" | sed 's/[^0-9]*//g' )



 if [ "$size_lincs" -lt "$size_pcts" ];
 then

  	size_train=$(  echo "$size_lincs * 0.8" | bc)
    size_test=$(  echo "$size_lincs * 0.2 + $size_train" | bc)
 else
    size_train=$(  echo "$size_pcts * 0.8" | bc)
    size_test=$(  echo "$size_pcts * 0.2 + $size_train" | bc)
 fi

size_train=$(echo ${size_train%%.*})
echo "Training size: $size_train"
size_test=$(echo ${size_test%%.*})
echo "Test size: $size_test"


# cut points
echo "PCT group: $pctGroup"
start_linc=$(  echo "2" | bc)
middle_linc=$(  echo "$size_train" | bc)
end_linc=$(  echo "$size_test" | bc)

start_pct=$(  echo "$size_test * $pctGroup + 2" | bc)
middle_pct=$(  echo "($size_test * $pctGroup) + 2 + $size_train" | bc)
end_pct=$(  echo "($size_test * $pctGroup) + 2 + $size_test" | bc)


echo "PCT  - start: $start_pct    middle: $middle_pct  end: $end_pct"
echo "LINC - start: $start_linc   middle: $middle_linc end: $end_linc"

resultFolder=$(  echo "result$pctGroup" )
echo "resultFolder: $resultFolder"


line_train_linc=$(  echo "$start_linc,""$middle_linc""p" )
line_test_linc=$(  echo "$middle_linc"",""$end_linc""p" )
line_train_pct=$(  echo "$start_pct,""$middle_pct""p" )
line_test_pct=$(  echo "$middle_pct"",""$end_pct""p" )
# echo "Group size: $tam2"

sed -n "1,1p"  "temp/linc_features.csv"> "temp/header.csv"
sed -n "$line_train_linc"  "temp/linc_features.csv"> "temp/lincs_train.csv"
sed -n "$line_test_linc" "temp/linc_features.csv" > "temp/lincs_test.csv"
sed -n "$line_train_pct"  "temp/pcts_features.csv" > "temp/pcts_train.csv"
sed -n "$line_test_pct" "temp/pcts_features.csv" > "temp/pcts_test.csv"
cat "temp/header.csv" "temp/lincs_train.csv" "temp/pcts_train.csv" > "temp/train.csv"
cat "temp/header.csv" "temp/lincs_test.csv" "temp/pcts_test.csv" > "temp/test.csv"


Rscript buildModel.r

mv "result" "$resultFolder"
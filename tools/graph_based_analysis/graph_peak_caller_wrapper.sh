source activate py36
pip install pyvg
exit 0

genome=$1
sample=$2
read_length=$3
fragment_length=$4
control=$5
out_fasta_file=$6
out_json_file=$7

echo "Running with sample $sample."

has_control="True"
if [ $control = "None" ]; then
    has_control="False"
    control=$sample
    echo "Will not use control (meaning sample will be used as control)"
else
    echo "Will use separate control track $control"
fi

# Get hardcoded genome size
if [ $genome = "drosophila_melanogaster" ]; then
    genome_size="98000000"
    chromosomes="chr3R,chr3L,chr2R,chr2L,chrX,chr4"
elif [ $genome = "arabidopsis_thaliana" ] ; then
    genome_size="135000000"
else
    >&2 echo "Error: Invalid genome."
    exit 1
fi

echo "Will use genome size $genome_size"

graph_dir="/home/ivargry/dev/galaxy/static/$genome/"
unique_reads=$(pcregrep -o1 '"sequence": "([ACGTNacgtn]{20,})"' $control | sort | uniq | wc -l)

echo "Found $unique_reads unique control reads."
echo "Has control $has_control"

echo "Splitting reads"
graph_peak_caller split_vg_json_reads_into_chromosomes $chromosomes $sample $graph_dir > log_splitting.txt 2>&1

if [ $sample != $control ]; then
    graph_peak_caller split_vg_json_reads_into_chromosomes $chromosomes $control $graph_dir
fi


if [[ $(wc -l <$control) -ge 2 ]]; then
    echo "Control is not empty"
else
    echo "Error: Control track is empty."
    exit 1
fi

sample_base_name=$(echo $sample | cut -f 1 -d '.')
control_base_name=$(echo $control | cut -f 1 -d '.')

for chromosome in $(echo $chromosomes | tr "," "\n")
do
    graph_peak_caller callpeaks_whole_genome $chromosome $graph_dir/ $graph_dir/ \
            $graph_dir/linear_map_ "${sample_base_name}_" "${control_base_name}_" "" $has_control \
            $fragment_length $read_length True $unique_reads $genome_size > log_before_pvalues_chr$chromosome.txt 2>&1 &
done

wait

for chromosome in $(echo $chromosomes | tr "," "\n")
do
	graph_peak_caller callpeaks_whole_genome_from_p_values \
                $chromosome $graph_dir/ "" $has_control $fragment_length $read_length > log_after_pvalues_chr$chromosome.txt 2>&1&
done

wait

graph_peak_caller concatenate_sequence_files $chromosomes $out_fasta_file > log_concatenating.txt 2>&1
cat *_max_paths.intervalcollection >> $out_json_file
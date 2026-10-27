# For unit testing
#counts_dir="/home/jeszyman/mpnst/analysis/cfdna-wgs/frag/counts"
#out_tsv="/home/jeszyman/mpnst/analysis/cfdna-wgs/frag/frag_counts.tsv"

# Define variables
counts_dir="${1}"
out_tsv="${2}"

# Remove the existing aggregate file if present
if [ -f $out_tsv ]; then rm $out_tsv; fi
#touch $out_tsv

# Make aggregate file
for file in ${counts_dir}/*;
do
    ## Add file name to each line
    #awk '{{print FILENAME (NF?"\t":"") $0}}' $file |
    #    # Modify file name to library id
    #    sed 's/^.*lib/lib/g' |
    #    sed 's/_.*_/\t/g' |
    #    # Cleanup "tmp"
    #    sed 's/.tmp//g' |
    #    # Send to output
    #    sed 's/\.bed//g' >> $out_tsv
    #
    # Extract the library name (e.g., healthy3 from the path)
    lib_name=$(basename "$file" | sed 's/_cnt_.*//g' | sed 's/\.tmp//g')

    #count length
    if [[ "$file" == *"short"* ]]; then
        len_class="short"
    elif [[ "$file" == *"long"* ]]; then
        len_class="long"
    else
        len_class="unknown"  # Fallback in case neither short nor long is found
    fi

    # Add file name to each line with library name and filter out any rows containing nan/NaN
    awk -v lib="$lib_name" -v len_class="$len_class" '
        tolower($0) !~ /nan/ {
            print lib"\t"len_class"\t"$0
        }' $file >> $out_tsv
   
    #awk -v lib="$lib_name" '{print lib"\t" $0}' $file |
    #    # Send to output
    #    sed 's/\.bed//g' >> $out_tsv
done

# Add a header
sed -i  '1 i\library	len_class	chr	start	end	gc	count' $out_tsv

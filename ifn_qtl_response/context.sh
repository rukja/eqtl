#!/bin/bash
#$ -cwd
#$ -o log_extract.txt
#$ -j y
#$ -l h_rt=2:00:00,h_data=2G
#$ -pe shared 3
#$ -m bea
#$ -M rithwiknarendra@ucla.edu

set -euo pipefail  # Improved error handling

pwd
. /u/local/Modules/default/init/modules.sh

module load qtltools/1.3.1
module load bcftools/1.11
module load htslib/1.12
module load R
module load vcftools

# Goals: The overall goal for this script is to create a new file containing context EQTLs. Each line should include information about a variant found in one group but not the other.
# For example, if a variant was only found in the VPA group, it would have its information followed by its group, which is VPA.
# Then, add the rsq and p value for this variant in the nominal file of the other group.


if [ "$#" -ne 1 ]; then
    echo "Usage: $0 condition"
    echo "Where condition contains a condition (ex:VPA)"
    exit 1
fi

var="$1"

base_dir=$(pwd)
final_output="${base_dir}/compare_unique_variants.txt"

echo "Starting processing. Base directory: ${base_dir}"

for chr in {1..22} X; do #For each chromosome (loop through 1-22 and X)
    dir="chr${chr}" #Set the directory based on the loop
    echo "Processing directory: ${dir}"

    if [ -d "$dir" ]; then # If the directory exists
        cd "$dir" # Enter the directory
        echo "Changed to directory: $(pwd)" # Print current path

        nom_file="filtered_nominals_${chr}.txt" # Establish control nominal file
        m_nom_file="${var}_filtered_nominals_${chr}.txt" # Establish treatment nominal file
        encyclopedia_nom="nominals_${chr}_1.txt" # Establish control "encyclopedia" file
        encyclopedia_nom_var="${var}_nominals_${chr}_1.txt" # Establish treatment "encyclopedia" file

        echo "Looking for files: ${nom_file} and ${m_nom_file}" # Next, check if nominal files exist

        if [[ -f "$nom_file" && -f "$m_nom_file" && -s "$nom_file" && -s "$m_nom_file" ]]; then # Check if the file exists and is a regular file (-f) and has a size greater than 0 (-s)
            echo "Found both files for chromosome ${chr}"

            header=$(head -n 1 "$nom_file") # Extract the header of the filtered control nominal file and set it as the header
            if [ "$chr" = "1" ]; then # Just for the first chromosome
                echo -e "$header group other_rsq other_pval" > "$final_output" # Print the header of the filtered control nominal file and other columns to be filled in to the final output. The -e flag enables interpretation of escape characters in each string (there are none, for now)
                echo "Created output file with header"
            fi

            echo "Processing chromosome ${chr}..."
            
            # Extract control-only variants in a single pass (replacing comm -23)
            awk 'BEGIN {FS= " "; OFS=" "}
                    # The OFS command within awk is used to set an output field seperater. It should be space delimited. We use begin because we want to do these actions before we # actually process the files
                 NR==FNR {gene_var[$1 FS $8 FS $3 FS $4] = $0; next}
                 # NR is the ordinal number of the current record from the start of input. The value is 0 inside a Begin function. FNR is the ordinal number of the current # # #record. Inside a Begin function, its value is 0. Put more elegantly, As FNR counts lines in each file and NR counts total lines, so (NR == FNR) is a test that # you are currently reading the first file named in the arguments.

                 # FS is an input field seperator (space by default). We create an associative array, storing the gene-variant pair ($1 and $8, separated by a space via FS) and # its corresponding line ($0). Simply put, the gene-varaint pair functions as the key and its line is the value. Next skips to the next line without finishing # the script. Because we set NR==FNR, this only happens for the first file.

                 # Now, we move onto the second file. For each line, if the gene and variant pair does not exist in our dictionary, print that line. This must be a treatment
                 # only variant

                 !($1 FS $8 FS $3 FS $4 in gene_var) {print $0, "'$var'", "placeholder1", "placeholder2" >> "'$final_output'"}' \
                "$nom_file" "$m_nom_file"
            
            # Extract treatment-only variants in a single pass (replacing comm -13)
            awk 'BEGIN {FS= " "; OFS=" "}
                 NR==FNR {gene_var[$1 FS $8 FS $3 FS $4] = 1; next}
                 !($1 FS $8 FS $3 FS $4 in gene_var) {print $0, "control", "placeholder1", "placeholder2" >> "'$final_output'"}' \
                "$m_nom_file" "$nom_file"
            
            # Process encyclopedia data for treatment variants
            if [ -f "$encyclopedia_nom_var" ]; then # First, check if the nominal encylopedia exists
                echo "Found variant encyclopedia file"
                awk 'BEGIN {FS= " "; OFS=" "}
                     NR==FNR && $19=="control" {control_variants[$1 FS $8] = NR; next} #Start by reading only the first file, which in this case is the unique variant output #  # file. Only if the 19th field is control (which is what NF does), build an associative array containing the gene-variant pair for
                     # control-only variants, setting the value as the current line number in the unique variant output file.
                     FNR>1 && ($1 FS $8 in control_variants) { # Move on to the second file and see if this gene-variant pair exists. If it does,
                         line_num = control_variants[$1 FS $8]; # Set its line number as the value for the array key

                         cmd = sprintf("sed -i '' \"%ds/placeholder1/%s/; %ds/placeholder2/%s/\" \"%s\"", line_num, $13, line_num, $12, ARGV[1]);

                            # Replace the - with 13 (p value) and 12 (rsq). Sed edits in-place
                         system(cmd);
                     }' \
                    "$final_output" "$encyclopedia_nom_var"
            fi
            
            # Process encyclopedia data for control variants if file exists
            if [ -f "$encyclopedia_nom" ]; then
                awk 'BEGIN {FS= " "; OFS=" "}
                     NR==FNR && $19=="'$var'" {var_variants[$1 FS $8 FS $3 FS $4] = NR; next}
                     FNR>1 && ($1 FS $8 FS $3 FS $4 in var_variants) {
                         line_num = var_variants[$1 FS $8 FS $3 FS $4];
                         cmd = sprintf("sed -i '' \"%ds/placeholder1/%s/; %ds/placeholder2/%s/\" \"%s\"", line_num, $13, line_num, $12, ARGV[1]);

                         system(cmd);
                     }' \
                    "$final_output" "$encyclopedia_nom"
            fi

            echo "Finished processing chromosome ${chr}"
        else
            echo "Missing one or both files in ${dir}"
            echo "nom_file exists: [$(test -f "$nom_file" && echo "yes" || echo "no")]"
            echo "m_nom_file exists: [$(test -f "$m_nom_file" && echo "yes" || echo "no")]"
        fi
        cd "$base_dir"
        echo "Returned to base directory: $(pwd)"
    else
        echo "Directory $dir not found"
    fi
done

echo "Processing complete! Final output is in: ${final_output}"
echo "Number of lines in output file: $(wc -l < "$final_output")"

#!/bin/bash

#!/bin/bash

input_vcf="snps_3_100.vcf"
output_vcf="snps_3_100_recode.vcf"

awk '
BEGIN { OFS="\t"; pos=1 }
{
    if ($1 ~ /^#/ ) {
        print $0  # Print header lines as is
    } else {
        $1 = 1      # Change chromosome ID to 1
        $2 = pos++  # Assign sequential position values
        print $0
    }
}' "$input_vcf" > "$output_vcf"

echo "Processing complete. Modified VCF saved as $output_vcf"


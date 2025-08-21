#!/bin/bash
set -e

echo "Updating system..."
sudo apt-get update

echo "Installing primary analysis tools..."
# Core NGS tools
sudo apt-get install -y \
    sra-toolkit \
    pigz \
    fastqc \
    multiqc \
    trimmomatic \
    hisat2 \
    samtools \
    subread

echo"Installation successful!"
cat <<EOF
To note: Qualimap must be installed via conda or manually
Reference:
- http://qualimap.conesalab.org/
- https://anaconda.org/bioconda/qualimap
EOF

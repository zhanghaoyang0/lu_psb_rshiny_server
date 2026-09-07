#!/bin/bash
cd ./temp

# Clean all files except the three specific CSV files
echo "Cleaning temp directory..."
echo "Preserving: genomic_del.csv, protein_del.csv, transcript_del.csv"

# Remove all files except the specified ones
find . -type f ! -name "genomic_del.csv" ! -name "protein_del.csv" ! -name "transcript_del.csv" -delete

echo "Cleanup completed!"
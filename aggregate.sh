#!/bin/bash

# Define the patterns for the text files
pattern1="ACC_*.txt"
pattern2="OMP_*.txt"
pattern3="cuda_*.txt"

# Define the output file
output_file="all_numbers_columns.txt"

# Clear the output file if it exists
> "$output_file"

# Write the headers to the output file
echo -e "$pattern1\t$pattern2\t$pattern3" >> "$output_file"

# Get the list of files for each pattern
files1=($pattern1)
files2=($pattern2)
files3=($pattern3)

# Get the maximum number of files to handle different file counts
max_files=${#files1[@]}
if [ ${#files2[@]} -gt $max_files ]; then max_files=${#files2[@]}; fi
if [ ${#files3[@]} -gt $max_files ]; then max_files=${#files3[@]}; fi

# Loop through each file index
for ((i=0; i<max_files; i++)); do
  # Read the first line of each file, if it exists
  number1=""
  number2=""
  number3=""
  if [ $i -lt ${#files1[@]} ]; then number1=$(head -n 1 "${files1[$i]}"); fi
  if [ $i -lt ${#files2[@]} ]; then number2=$(head -n 1 "${files2[$i]}"); fi
  if [ $i -lt ${#files3[@]} ]; then number3=$(head -n 1 "${files3[$i]}"); fi
  
  # Append the numbers to the output file
  echo -e "${number1}\t${number2}\t${number3}" >> "$output_file"
done

echo "All numbers have been collected into $output_file in three columns"


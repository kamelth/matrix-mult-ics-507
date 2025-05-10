#!/usr/bin/env bash
set -euo pipefail

# Experiment parameters
methods=(Sequential SequentialP StraightDivAndConq StraightDivAndConqP StrassenDivAndConq)
ns=(2048 4096 8192 16384)
threads=(1 2 4 8)
bases=(16 32 64 128 256 512 1024 2048)

# Calculate total number of tasks
total=0
for method in "${methods[@]}"; do
  if [[ "$method" == "Sequential" || "$method" == "SequentialP" ]]; then
    total=$((total + ${#ns[@]} * ${#threads[@]}))
  else
    total=$((total + ${#ns[@]} * ${#threads[@]} * ${#bases[@]}))
  fi
done

# Output CSV
out=results.csv
echo "method,n,threads,base_thresh,time_s" > "$out"

# Initialize counter
counter=0

echo "Starting matrix multiplication experiments..."
for method in "${methods[@]}"; do
  for n in "${ns[@]}"; do
    infile="./inputs/input_${n}.txt"
    name="input_${n}"
    for t in "${threads[@]}"; do
      if [[ "$method" == "Sequential" || "$method" == "SequentialP" ]]; then
        ./matmul -i "$infile" -m "$method" -t "$t"
        info_file="outputs/${name}_${n}_info_${method}.txt"
        time_s=$(<"$info_file")
        echo "${method},${n},${t},,${time_s}" >> "$out"

        # Update progress with current method
        counter=$((counter + 1))
        percent=$(( counter * 100 / total ))
        echo -ne "\r[${method}] -n ${n} -t ${t} - Progress: ${percent}% (${counter}/${total})"
      else
        for b in "${bases[@]}"; do
          ./matmul -i "$infile" -m "$method" -b "$b" -t "$t"
          info_file="outputs/${name}_${n}_info_${method}.txt"
          time_s=$(<"$info_file")
          echo "${method},${n},${t},${b},${time_s}" >> "$out"

          # Update progress with current method
          counter=$((counter + 1))
          percent=$(( counter * 100 / total ))
          echo -ne "\r[${method}] -n ${n} -t ${t} -b ${b} Progress: ${percent}% (${counter}/${total})"
        done
      fi
    done
  done

done

echo -e "\nAll experiments done. Results saved to $out"

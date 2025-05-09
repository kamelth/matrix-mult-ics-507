# Matrix multiplication ICS-507

# Instructions to run the code

### 1. Generate input file
`python generate_input.py <n> ./inputs/<output_file_name>.txt`

### 2. Compile c code
`make all`

### 3. Select and run a single method
```bash
# single‐threaded Strassen on input1.txt
./matmul -i inputs/input1.txt -m StrassenDivAndConq

# 8 threads, cutoff=32 for parallel D&C
./matmul -i inputs/input1.txt -m StraightDivAndConqP -b 32 -t 8
```

### Possible parameters

| Flag                    | What it controls                                                   | Example values                                                                           |
|-------------------------|---------------------------------------------------------------------|------------------------------------------------------------------------------------------|
| `-i`                    | Input file                                            | `input1.txt`, `input2.txt`, …                                                                         |
| `-m`                    | Which algorithm to run                                              | `Sequential`, `SequentialP`, `StraightDivAndConq`, `StraightDivAndConqP`, `StrassenDivAndConq` |
| `-b`                    | Base-case size for divide & conquer (sequential or parallel)        | Any power-of-2 cutoff (e.g. 32, 64, 128)                                                  |
| `-t`                    | Number of OpenMP threads                                            | 1, 2, 4, 8, 16, …                                                                         |
| `OMP_NUM_THREADS` (env) | Fallback thread count when `-t` is not provided                     | 1, 2, 4, 8, 16, …                                                                         |

# README

## Overview

This repository contains the source code for our ongoing research on efficient quasi-clique detection in dynamic graphs.
Our approach focuses on handling graph updates (edge additions and deletions) efficiently, aiming to maintain and update large quasi-clique structures without recomputing from scratch after each change.
The code is written in C++, using only standard libraries, and can be compiled in any environment.

## Dataset

All datasets used in our experiments are mainly derived from **SNAP** and **KONECT**, with randomized edge additions and deletions to simulate dynamic updates.

There are five types of datasets:

- **DataSet_Only_Add**  
  Based on static datasets from SNAP, with only edge additions applied.

- **DataSet_Only_Delete**  
  Based on static datasets from SNAP, with only edge deletions applied.

- **DataSet_Random_Operation**  
  Based on static datasets from SNAP, with randomized edge additions and deletions.

- **DataSet_Temporal_Add_Deletion**  
  Derived from temporal datasets in KONECT, with randomized edge deletion times.

- **DataSet_Temporal_Incremental**  
  Derived from temporal datasets in KONECT, using their natural incremental evolution.
  
## How to run the code

To compile the program, put all source files in the same directory and simply run:

```bash
make
```

Then you can execute "DynamicQC". The format of the input parameter is

```bash
./DynamicQC -f filepath -gamma para1 -b para2 -k para3 -l para4 -alpha para5 -Batch para6 -B para7 -Q para8
```

For example, `./DynamicQC -f DataSet_Random_Operation/facebook.txt -gamma 0.9 -b 0.6 -k 8 -l 8 -alpha 0.9 -Batch 5000 -Q 100000` is running `DynamicQC` on facebook dataset (with randomized edge updates) using the following parameters: gamma = 0.9, b = 0.6, k = 8, l = 8, alpha = 0.9, Batch = 5000, B = 5 and Q = 100000.

If you want to execute `DynamicQC` with BottomK method, the fomat of the imput parameter is,

```bash
./DynamicQC -f filepath -gamma para1 -b para2 -k para3 -l para4 -alpha para5 -Batch para6 -B para7 -Q para8 -use_bottom_k
```

### Running the baseline method (NSF)

We also include the implementation of the baseline method NSF for comparison.

To compile NSF, use the following command:

```bash
g++ NSF.cpp NSF.h -o NSF
```

Then you can execute `NSF`.
For the `NSF-fn` variant, the input format is:

```bash
./NSF input-filename gamma b k
```

For the `NSF-ns` variant, the input format is:

```bash
./NSF input-filename gamma b
```

All methods use the same input data format to ensure a fair comparison.

Unless otherwise specified, all experiments use identical parameter settings and update streams across different methods.

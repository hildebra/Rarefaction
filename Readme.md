# Rarefaction tool kit - rtk
A rarefaction software written in C++11 to rarefy large high count datasets quickly and return diversity measures.

# Installation
To use rtk you can download the binary files under https://github.com/hildebra/Rarefaction/releases or compile from source.

## Compile from source
To build this software you will need to have a compiler for C++11 on your system. On a GNU/Linux system you usuallu have to install developer tools to do that. For Ubuntu this is explained here: https://help.ubuntu.com/community/InstallingCompilers


rtk was tested to compile successfully on windows (version? with compiler XY), GNU/Linux (g++ v. 4.8.5 and v. 6.1.1) and on Mac OS 10.11.2 (Apple LLVM version 7.0.0 (clang-700.0.72)).


**Compile in UNIX**
```bash
git clone https://github.com/hildebra/Rarefaction
cd Rarefaction/rtk
make
```

# Usage
Two modes for rarefaction of a count table are available

```bash
rtk  -i <input.csv> -m <mode> -o <output> [options]
```

### Options:
```
-i      path to an .csv file to rarefy
-o      path to a output directory
-m      mode can be either swap or memory.
        Swap mode creates temporary files but uses less memory.
        The speed of both modes is comparable.
-d      Depth to rarefy to. Default is 0.95 times the minimal column sum.
-r      Number of times to create diversity measures. Default is 10.
-w      Number of rarefied tables to write.
-t      Number of threads to use. Default: 1
-ns     If set, no temporary files will be used when writing rarefaction tables to disk (no swap).

```

### Output files:

**median_alpha_diversity.tsv**

This file contains the median diversity measures for all Samples in a tab separated format.

**[samplename]_richness|eveness|...|.tsv**

Each diversity measures is exported as a table containing all repeats for all sample.

**rarefied_to_X_n_Y.tsv**

If `NoOfMatrices > 0` each rarefied matrix will be saved in the output directory under this file. The structure of all files is the same and similar to the input file.

**sums.txt**

This file contains the column sums of all samples. It can be used to estimate well suited rarefaction depth.

### Temporary files
If the mode `memory` is used, temporary files will be produced to reduce RAM usage. Thus the input matrix will be first split into its columns and each column will be written into a single file. Those file will then be loaded again and deleted after the software is finished using them.

Temporary files will also be created if `-w > 0`. In this case the vectors of the rarefied tables will be stored on disk as binary before merging them to tables. Thi can be prevented by using the `-ns` flag.

In both cases RAM usage is drastically reduced and the load on the local drive is substantially higher.



## Input data format
Input data for rtk should be a count table in a .tsv or .csv format.
Row and column names must be provided and be unique.

**Example file:**

|       | Sample a | Sample b | Sample c | Sample d |
|-------|----------|----------|----------|----------|
| OTU 1 | 0        | 12       | 4        | 80       |
| OTU 2 | 5        | 30       | 0        | 10       |
| OTU 3 | 110       | 0        | 1       | 0        |
| OTU 4 | 43       | 253        | 15       | 30        |
| OTU 5 | 0       | 0        | 15       | 0        |
| OTU ... | ...       | ...        | ...       | ...        |
| OTU ... | ...       | ...        | ...       | ...        |
| OTU n | 25       | 12        | 3       | 0        |

Rarefaction is always performed on the columns of the dataset. If you want to rarefy on the rows please consider transposing your input data ahead of rarefaction.

### Transposing input data
On a UNIX system use AWK to transpose a `.csv` table:
http://stackoverflow.com/questions/1729824/transpose-a-file-in-bash





# Example
A minimal working example of a rarefaction is shown here. This example should run on any UNIX system.
```
#!/bin/sh
FILE="example.input.csv"
touch $FILE
echo -e "OUT    \tSample 1\tSample 2\tSample 3"       >> $FILE
echo -e "OTU 1\t  232      \t  10       \t  0"        >> $FILE
echo -e "OTU 2\t  0        \t  57       \t  22"       >> $FILE
echo -e "OTU 3\t  17       \t  0        \t  45"       >> $FILE
echo -e "OTU 4\t  5        \t  83       \t  0"        >> $FILE

./rtk -i $FILE -o test. -mode memory
ls -lh test.*
```


# Copyright
rtk is licensed under the GPLv2. See notice and license file for more information.

Copyright (c) 2016 by Falk Hildebrand and Paul Saary

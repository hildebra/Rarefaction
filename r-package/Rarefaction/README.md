# Rarefaction - R module
This R module is a wrapper for the rarefaction programm.

## Dependencies
This package depends on the Rpackage RCPP, which enables using Cpp-Code easily in R 
packages.

## How to install
To build and install this package from source, just run these commands in a directory 
above the dir containing the source files.
```
R CMD build Rarefaction
R CMD install Rarefaction_*.tar.gz
```

## Running the package
```R
require('Rarefaction')
path <- '/path7to/a/csv.file
rare(
```

## State of the project
The project is currently developed in the devlopment branch.
No stable version is currently provided.

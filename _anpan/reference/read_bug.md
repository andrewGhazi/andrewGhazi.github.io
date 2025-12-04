# Read a bug's genefamily file

Read a bug's genefamily file

## Usage

``` r
read_bug(bug_file, meta = NULL, remove_pattern = "_Abundance-RPKs")
```

## Arguments

- bug_file:

  path to a bug's genefamily file

- meta:

  a data frame of metadata

- remove_pattern:

  pattern to remove from the column names of the genefamily file

## Value

The genefamily file of the bug as a data.table in TALL format

## Details

The input bug_file needs to be readable by data.table::fread()

If metadata is provided, the genefamily file is subset to only those
samples present in the sample_id column of the metadata.

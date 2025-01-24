1. Prepare fasta files for query genome and target genome in data directory
2. Change the sequence name (the line starting with > in the fasta file) as you want
3. write chroms.txt file which match the chromosome name of query name and target name. For exmpale,
```
Chr1,chr1
Chr2,chr2
Chr3,chr3
Chr4,chr4
Chr5,chr5
```
The chromosome name should be listed in the order of query then target, and the two should be separated by comma.

4. edit config.yaml file accordingly.

Now you're all set.

`snakemake -p --cores <number of threads>`

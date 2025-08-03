close all; clear all;

file_gwas='GWAS_mtcorr_pgen.phenotype.glm.linear';

t=readtable(file_gwas,'FileType','delimitedtext');

etc_gwas_manhattan(t);
# **Matrix genartion**

The genomic folder consists of: </p>
-          **input** folder with test input files,
-          **Reference** folder with all references,
-          **output** folder,
-          **scripts** folder.

The code allow you to generate the matrix for the prediction model algorithm. </p>

The scripts are 5:
- 1.CNV_manage.R -> this code consists of two parts.
The first is meant to create a matrix with the GISTIC peak alterations that we found on the training set (N=1933). At beginning we implemented a distinction between focal and large events (it is present in this code and it is used on a different analysis in the paper) but at the end we don't consider this difference (see paper and Supplementary Data 1), in fact this information disappears in the matrix (using the next codes). </p> 
One important thing to check is the copy-number values. We use rounded integers, so gains are =3, amp =4. If you have non-rounded value probably you will have to change these cut-off (lines: 119:122 and 170:176).
 
The second part is meant to create a matrix for the hyperdiploid chromosomes. We considered hyperdiploid when 60% of the arms was at least with a gain CN. We did not consider the duplicated genomes.
 
Brief summary of each R file:
- 2.Cytogenetic_abnormalies.R -> with this code you can create the HRD, other GISTIC peak, 1q peak and translocation matrix.
- 3.Mutated_genes.R -> you can create the ONCOGENES, Tumor Suppressor Genes (TSGs), OTHER genes and CHRX (genes in chrX, we consider only - mutations and no copy-numbers) matrix.
- 4.MutationalSignature.R -> to identify APOBEC using mmsig.
- 5.Matrix_generation.R -> to assemble the final matrix using the intermediate files previously created.
 
When you run all these codes, the intermediate outputs will be created in the **output** folder, while the final matrix will be generated in the "genomic" folder.

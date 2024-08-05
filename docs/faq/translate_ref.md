# Finding a `translate` and `align` reference file

In order to correctly align and annotate the genes in your sequences, augur uses an annotated reference sequence in 'Genbank' format (if Fasta input) or 'GFF' format (if VCF input).

You should find an appropriate reference file that corresponds to your sequence data. For Fasta input, you'll then use this sequence for the `align` step and the `translate` step. (It's important you use it for both so that the numbering is consistent.)

### Finding the right sequence

You can look for a reference sequence on Genbank. It should:
* Include all the genes you wish to annotate (as individual annotations)
* Be approximately the same length as your sequence (Fasta-input)
* Not be too divergent from your sequence (or it may have trouble aligning) (Fasta-input)
* Use the same numbering as the reference your VCF is mapped to (VCF-input)

For example, in our Zika build we use the reference sequence [PF13/251013-18](https://www.ncbi.nlm.nih.gov/nuccore/KX369547). Since we have full-genome Zika samples, it was important we choose a reference that is also full-genome in length. If you have partial-genome sequences, you should look for a reference that covers the same region you have.

### Fasta

You can download a reference sequence you find on Genbank by using the 'send to' link in the top-right corner and then selecting 'Complete Record', 'File', and 'GenBank' as the format.

You'll then need to modify the Genbank file. We recommend not translating things like 'polyproteins' which contain many genes, and instead translating individual genes. Remove any genes/CDS which you do not want to translate.

It is important you leave the `source` information.

Augur will only translate genes which have 'CDS' as the feature name, and have an attribute called 'gene' or 'locus_tag'. You may have to modify your Genbank file so that it's in the correct format. Here's an example:

```
     CDS             961..2472
                     /product="envelope protein"
                     /gene="ENV"
```

Compare the [original Zika reference on Genbank](https://www.ncbi.nlm.nih.gov/nuccore/KX369547) to the [one used on Nextstrain](https://github.com/nextstrain/zika/blob/-/phylogenetic/defaults/zika_reference.gb). Notice that the genes are designated by `CDS` instead of `mat_peptide` and have an entry for `/gene=` as well as `/product=`.

### VCF

You can also find an appropriate GFF annotation reference on GenBank. Be sure to pick one that is very close to the strain you are using - especially if there might be variability in the genes present! If the positions in the GFF file do not match the positions in your VCF file, it will not work.

To download a file from GenBank, use the 'send to' button in the top-right corner, then select 'Complete Record', 'File', and 'GFF3' for the format.

You will have to modify the first column of the GFF so that it matches the `CHROM` (first column) of your VCF file.

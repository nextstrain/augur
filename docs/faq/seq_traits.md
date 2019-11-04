# Inferring Sequence Traits (like Drug Resistance)

The augur function `sequence-traits` can identify any trait associated with particular nucleotide or amino-acid mutations, but it's often used to identify drug resistance mutations (DRMs).

To tell augur which sites confer what trait (or drug resistance), you'll need to pass a file detailing these sites. 

<span style="color:red">TODO: Doesn't have to contain all 5 columns! Update to reflect!</span>

The file should contains five columns: GENE, SITE, ALT, DISPLAY_NAME, and FEATURE. DISPLAY_NAME can be blank.

### Amino Acid Sites

For example, for drug resistance in TB, we list the gene, the AA position in the gene, the AA mutation that confers resistance (you can list a site multiple times if multiple bases give resistance), and the name of the drug this mutation gives resistance to:

```
GENE    SITE    ALT DISPLAY_NAME    FEATURE
gyrB    461 N       Fluoroquinolones
gyrB    499 D       Fluoroquinolones
rpoB    432 E       Rifampicin
rpoB    432 K       Rifampicin
```

We can leave DISPLAY_NAME blank, as auspice will by default display the gene, site, and original and alternative base.

### Nucleotide Sites

For mutations outside of protein-coding genes, we can specify their position using nucleotides:

```
GENE    SITE    ALT DISPLAY_NAME    FEATURE
nuc 1472749 A   rrs: C904A  Streptomycin
nuc 1473246 G   rrs: A1401G Amikacin Capreomycin Kanamycin
nuc 1673423 T   fabG1: G-17T    Isoniazid Ethionamide
nuc 1673425 T   fabG1: C-15T    Isoniazid Ethionamide
```

In the TB literature, these mutations are still referred to by their position within non-protein-coding genes (`rrs`) or location near genes (`-17 fabG1`), not their nucleotide location. We can ensure auspice displays the more useful common nomenclature by giving entries for the DISPLAY_NAME column.

### Options

`sequence-traits` will return a value for each "feature" - for example, all the mutations on the tree that lead to resistance to Streptomycin. It will also generate a count either of the total number of "features" each node has (ex: the total number of drugs a sequence is resistant to), or the total number or mutations specified in the file each node has (ex: the total number of DRMs a sequence has, even if some are for the same drug). 

You can specify a name for this count using the `--label` argument (ex: "Drug_Resistance"). The `--count` argument value specifies whether to count the number of traits (ex: drugs resistant to) (use `traits`) or number of overall mutations (use `mutations`).


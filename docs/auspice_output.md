
## Format of JSON outputs for Auspice

#### tree JSON
Augur outputs the tree as a JSON file modeled after the format used by d3.
```
{
    strain: 'name of virus or internal node',
    clade:  1,       //clade identifier, typically an integer
    xvalue: 0.052,   //tree layout: divergence from root
    yvalue: 251,     //tree layout: vertical position in tree
    tvalue: 6.25     //tree layout: years since root, num_date attribute can be used instead
    muts: [ ('A135G'), () ...],  //nucleotide mutations mapped to branch
    aa_muts: {       //amino acid mutations mapped to branch
        'protein1':[ ('F155Y'), () ...]},
        'protein2':[]
        },
    attr: {
        region: "africa",   //top level geo
        country: "nigeria",
        city: "lagos",
        num_date: 2015.14,  // numerical date
        date: 2015-02-13    //YYYY-MM-DD date string
        div: 0.052,         //divergence from root
        // any other attributes that auspice needs to know
    }
    children:[
        {
            strain: ...
        },
        {
            ...
        },
    ]
}
```

#### sequence JSON
The input sequences are stored in a compressed JSON format as follows:
```
{
    root: { //the sequence inferred for the root node
        nuc: "ACGAGTGATG...",
        protein1: "KCYTWD...",
        protein2: "ASRTRTY...",
    },
    1: {   // sequences for each node in the tree (incl the root)
           // these can either be full sequences or differences relative to root
        nuc: { 135:"G", 225:"A"},
        protein1: {155: "Y"},
        ...
    },
    2: {
        ...
    },
    ...
}
```

#### frequency JSON
```
{
    pivots: [2013, 2013.25, 2013.5 ...],     //interpolation pivots for frequency trajectories
    // identifier format: region_category:position
    "global_protein1:135Y"  : [0.0, 0.1, ...],    //frequencies at the pivots points
    "global_nuc:135G"       : [0.0, 0.1, ...],    //nucleotide mutation frequencies
    "global_clade:2"        : [0.0, 0.1, ...],   //clade frequencies
}
```

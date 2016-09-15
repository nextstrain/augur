## Introduction

The nextstrain project is an attempt to make flexible informatic pipelines and visualization tools to track ongoing pathogen evolution as sequence data emerges. The nextstrain project derives from [nextflu](https://github.com/blab/nextflu), which was specific to influenza evolution.

nextstrain is comprised of three components:

* [fauna](https://github.com/nextstrain/fauna): database and IO scripts for sequence and serological data
* [augur](https://github.com/nextstrain/augur): informatic pipelines to conduct inferences from raw data
* [auspice](https://github.com/nextstrain/auspice): web app to visualize resulting inferences

## Augur

*Definition: One held to foretell events by omens.*

Augur is the informatic processing pipeline to track evolution from sequence and serological data.  It aims to

* subsamples, cleans and align sequences
* build a phylogenetic tree from this data
* infer ancestral states
* infer timing of phylogenetic branching events
* infer mutation and clade frequency trajectories through time
* export a JSON bundle for visualization

### Format of JSON outputs

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
    global_protein1:135Y: [0.0, 0.1, ...],    //frequencies at the pivots points
    global_nuc:135G: [0.0, 0.1, ...],         //nucleotide mutation frequencies
    global_clade:2:  [0.0, 0.1, ...],         //clade frequencies
}
```

## License and copyright

Copyright 2014-2016 Trevor Bedford and Richard Neher.

Source code to nextstrain is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). nextstrain is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.

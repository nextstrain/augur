## AUGUR
Augur consists of 2 stages: prepare & process. You prepare the data (from fauna, most commonly) into a series of JSONs - one per segment. Each of those can then be processed independently.


## Docs:
* [Prepare](prepare.md)
* [Process](process.md)
* [Format of (auspice) output JSONs](auspice_output.md)

### Status of branch `rejig`:

| Virus         | Status           |
| ------------- | ------------- |
| Seasonal flu      | partially working - see flu/README.md |
| Zika    | working      |
| Ebola | working      |
| Dengue | partially working (no freqs, no titers)      |
| H7N9 | working      |

Functionality missing everywhere:
* Reference sequences are not included in the analysis
* Restore ability limited
* No frequencies
* No titers

## How to run
* see the README.md files in the respective pathgen's folder

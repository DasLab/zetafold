
#zetafold/analysis
## What this is
Some analysis scripts to compare zetafold's accuracy across its training and test sets with different parameters, and also with other packages.

## Output
Zetafold or the other packages are run from command line, with results stored in subdirectories of your run directories, with the *base pair probability file* (e.g., `bpp.txt.gz`) being particularly important.

Then you'll create an analysis table with:

 • **Ensemble defect** as defined in NUPACK. Average of the 'defect' over all nucleotides. If the nucleotide is supposed to be unpaired, consider its probability of being paired; if the nucleotide is supposed to be paired to a particular partner, pairing to other partners or being unpaired is a defect.  

 • **BPP defect** Just scores the defects as nucleotides that are supposed to be base paired. The idea is that for many databases of 'native structures' (e.g., RGAM), paired regions are inferred from phylogenetic covariance, unpaired regions are actually unknown and might involve base pairs. 
 
##Installation
Make sure to include in your path this directory (analysis), e.g., with a command like this in your `.bashrc`:

```
PATH=$PATH:$HOME/src/zetafold/analysis/
```

and add the zetafold repository to your python path, e.g., something like

```
export PYTHONPATH=$PYTHONPATH:$HOME/src/zetafold/
```

also in your `.bashrc`. Be sure to update the paths above to point to your actual zetafold locations!

## Example

To run packages and accumulate output in subdirectories:

```
run_packages.py  --data favorites --packages zetafold_v0.171 -j4
run_packages.py  --data favorites --packages zetafold_v0.18 -j4
```

(the '-j4' means use 4 cores)

then run:

```
create_analysis_table.py  --data favorites --packages zetafold* 
```

to get an analysis table that should look like:

```
Ensemble defect
                           RNA zetafold_v0.171  zetafold_v0.18
                          tRNA          0.0277          0.0235
                          P4P6          0.2791          0.2444
        add_riboswitch_aptamer          0.0382          0.0270
                       5S_rRNA          0.6104          0.6205
     cdiGMP_riboswitch_aptamer          0.3287          0.3340
 FN_glycine_riboswitch_aptamer          0.4045          0.3766
     AdoCBL_riboswitch_aptamer          0.2064          0.2079
        FMN_riboswitch_aptamer          0.3674          0.3375
             HCV_IRES_domainII          0.4740          0.4307
                       OVERALL          0.3209          0.3052
                           RNA zetafold_v0.171  zetafold_v0.18

BPP defect
                           RNA zetafold_v0.171  zetafold_v0.18
                          tRNA          0.0353          0.0292
                          P4P6          0.2593          0.2075
        add_riboswitch_aptamer          0.0478          0.0282
                       5S_rRNA          0.7106          0.7135
     cdiGMP_riboswitch_aptamer          0.3919          0.3980
 FN_glycine_riboswitch_aptamer          0.5194          0.4729
     AdoCBL_riboswitch_aptamer          0.2545          0.2548
        FMN_riboswitch_aptamer          0.4918          0.4352
             HCV_IRES_domainII          0.6457          0.6056
                       OVERALL          0.3806          0.3565
                           RNA zetafold_v0.171  zetafold_v0.18           
```

## Other packages

Also can supply `contrafold` and `contrafold-nc` (contrafold with noncomplementary parameter set, i.e. allowing noncanonicals) 
as "--packages" for 
`run_packages.py`, e.g.: 

```
run_packages.py  --data favorites --packages contrafold contrafold-nc
```

You'll need to make sure the contrafold exectuable directory (usually `contrafold/src/`) is in your path, i.e. there should be a line in your .bashrc like this:

```
PATH=$PATH:$HOME/src/contrafold/src/
```

TODO: Vienna, RNAstructure, RNAsoft, etc.


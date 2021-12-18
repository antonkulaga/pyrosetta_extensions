# pyrosetta_extensions #

pyrosetta is a very powerful library with super-ugly API. 
This repository provides helper methods/classes to make it more usable

## setting things up ##

Pyrosetta has dual licensing policy, so check that your use-case is non-profit, see more info at http://pyrosetta.org/downloads
Note: conda version of pyrosetta does not include some code, so it is recommended to download wheel instead

```bash
micromamba env create -f environment.yaml
micromamba activate pyrosetta_extensions
```
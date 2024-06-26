This folder contains the scripts that build the curated datasets in `../curated-data` from the original datasets in `../raw-data`. 


- `build-curated-datasets.R` will invoke the core scripts in `GBI` and `TLI` folders, building the curated datasets
- `functions.R` contains the functions shared by the scripts in this repository
- `lang-metadata.csv` contains the taxonomic and geographic metadata on languages, largely from Glottolog. Missing coordinates and missing assignments to Glottolog macro areas were manually added. Assignments to AUTOTYP continents and areas were done manually.
- `GBI` contains the scripts and configuration data required to build the curated `GBI` datasets
- `TLI` contains the scripts and configuration data required to build the curated  `GBI` datasets

# Hantavirus Mutations Tracker

Purpose: This Hantavirus Mutations Tracker is a Python-based bioinformatics pipeline that retrieves global Hantavirus genomic data, aligns variants using remote NCBI BLAST, and locates mutational hotspots across the S-segment genome.

## Description

Hantaviruses are RNA viruses that infect rodents and are occassionally transmitted to humans. They are the major cause of two global diseases: hantavirus pulmonary syndrome (HPS), or hantavirus cardiopulmonary syndrome (HCPS), and hemorrhagic fever with renal syndrome (HFRS). 

The genome of the Hantavirus is comprised of three negative-sense, single-stranded RNAs which are highly conserved. 
* The large segment (L segment) encodes the RNA-dependent RNA polymerase (RdRp) responsible for transcription and replication of the viral genome. 
* The medium segment (M segment) encodes two surface glycoproteins (Gn and Gc). Gn and Gc are the primary surface glycoproteins that play a critical role in virus cell entry and virus assembly.
* The small segment (S segment) encodes the nucleoplasmid (N) protein. The N protein is involved in genome packaging, intracellular protein transport, and other crucial processes during hantavirus infection.

This project isolates and analyzes the small (S) segment.

## Getting Started

### Dependencies
1. Ensure you have [Python 3.8+](https://www.python.org/downloads/) installed.
2. Install [Biopython](https://biopython.org/)
3. Install [pandas](https://pandas.pydata.org/)
4. Install [matplotlib](https://matplotlib.org/)
5. Install [seaborn](https://seaborn.pydata.org/)

### Installing
1. Clone the repository.
``` bash
git clone [https://github.com/renee-do/Hantavirus-mutations-tracker.git](https://github.com/renee-do/Hantavirus-mutations-tracker.git)
cd hantavirus-mutations-tracker
```

## Authors

Renee Do
[LinkedIn](https://www.linkedin.com/in/renee-do)
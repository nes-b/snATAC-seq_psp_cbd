# Single-nucleus ATAC-sequencing in PSP and CBD brains


## Background

This repository contains scripts and metadata of our human brain tissue study: 

**Single Nucleus Accessible Chromatin Profiling Highlights Distinct Astrocyte Signatures in Progressive Supranuclear Palsy and Corticobasal Degeneration**. 
*https://doi.org/10.1007/s00401-022-02483-8*


### Prerequisites

- Computing environment with >= 200 GB RAM and >= 12 CPU cores for the pre-processing (01_cellranger_snappre.Rmd) 
- Computing environment with >= 16 GB RAM (+ extra swap-memory) and >= 12 CPU cores for all other steps (recommended, Linux-based distro recommended)  
- RStudio Server or Rstudio Desktop running with R3.6 or R4.0 (for *pathfindR*)
- ImageJ for running image pre-processing macros

Please see https://www.r-project.org/ and https://imagej.nih.gov/ij/index.html

- Raw data (bam files) can be accessed via the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/view/PRJEB54978)

### Code Execution

The code can be run from within R and ImageJ/Macros.


## Contributing

Dr. Viktoria C Ruf, Katrin Pratsch, Dr. Sigrun Roeber, Jeannine Widmann, Janina Mielke, 
Dr. Dr. Mario M Dorostkar, Dr. Otto Windl, Dr. Thomas Arzberger, 
Prof. Dr. Jochen Herms*, [Dr. Felix L Strübing](https://github.com/fstrueb)*

* contributed equally

Center for Neuropathology, \
German Center for Neurodegenerative Diseases, Translational Research \
University Hospital Munich \
Ludwig-Maximilians-University \
Germany


## Maintenance Git-Hub Repository

* **Nils Briel** - [nes-b](https://github.com/nes-b)

See also the list of [contributors](https://github.com/nes-b/snATACseq_psp_cbd/blob/master/contributors.txt) who participated in this project.


## Citation

If you use parts of this workflow, please cite:
- the original article: Briel, N., Ruf, V.C., Pratsch, K. et al. Single-nucleus chromatin accessibility profiling highlights distinct astrocyte signatures in progressive supranuclear palsy and corticobasal degeneration. Acta Neuropathol (2022). https://doi.org/10.1007/s00401-022-02483-8

- This repository: https://github.com/nes-b/scATACseq_psp_cbd 


## Ethics statement

The available resources stem from deceased individuals, who had confirmed to the use of their brain tissue for biomedical research.
All data is anonymised, so that access to the actual identifiers is restricted to the contributors themselves. \
This project had been accepted by the Ethic Board of the Ludwig-Maximilians-University, Munich.


## License

This project is licensed under the GNU License (version 3.0) - see the [LICENSE](LICENSE) file for details


## Acknowledgments

* Thanks to the members of the Center for Neuropathology, the Neurobiobank Munich and German Center for Neurodegenerative Diseases for their valuable conceptional and technical input.
* We would like to express special thanks to J.M. Luque, K. Pratsch, J. Widmann, and K. Ochs for fruitful discussions. 
* Furthermore, our thanks are addressed to Prof. M. Kerschensteiner, Dr. E. Beltran and their lab members for their structural support.
* Finally, we are greatful for all patients and their family members to enable this research.

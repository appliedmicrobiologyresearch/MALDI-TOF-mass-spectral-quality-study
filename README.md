## Factors associated with MALDI-TOF mass spectral quality of species identification in clinical routine diagnostics

MALDI-TOF MS is the most commonly used tool for microbial species identification. One of the challenges in microbial species identification by MALDI-TOF MS is mass spectral quality, which is currently vaguely defined. 
In this project we identify MALDI-TOF mass spectral quality features which are associated with a correct species identification by three different databases. We use these to assess the impact of varying measurement conditions on mass spectral quality.

All the raw- and metadata required to run these scripts are available at the Open Science Foundation (https://osf.io/ksz7r/). 

The scripts 'eval_ascii_axima.R', 'eval_ascii_microflex.R' and 'eval_ascii_microflex_calibration.R' read out mass spectral features for spectra acquired either on an Axima Confidence or a microflex Biotyper MALDI-TOF MS device. The input files for these scripts are ASCII files including the m/z and intensity values of picked peaks.

The scripts 'readout_VitekMS_ID.R', 'readout_microflex_Biotyper_output.R' and 'summarise_marker_based.R' summarise the species identification yielded by comparison to the VitekMS, the MALDI Biotyper and the PAPMID database, respectively. 

The script 'spectra_comparison.R' summarises the spectra evaluation and species identification and produces a series of plots, comparing the mass spectra acquired under varying sample preparation protocols and from multiple phylogenetic groups. 
The script 'plot_which_endpoints.R' assesses which mass spectral features are associated with a correct species identification and therefor good proxies for mass spectral quality. 


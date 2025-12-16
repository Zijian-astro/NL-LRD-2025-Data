# Stacking the X-Ray signal of Little Red Dots (Yue+24)

This repository holds the code and data used in Yue+24 (link). The repository is organized as follows:

code/: includes all the code used in this work. The most important files are the following two:

- code/read_cstack_results.py: provides core functions for reading cstack outputs, computing and saving probability distributions of count rates, and producing the final table that contain all information for follow-up analysis. Additional deriviation can be found in (file).

- code/makeplots.py: produces figure 1 and 2 in the paper.

The following files contain utilities used for the analysis:

- code/nh_selenium.py: a script to obtain the Milky Way column density from the web.

- code/pimms_selenium.py: a script to obtain the conversion factor between count rates and physical fluxes from the web.

- code/convfactor.py: a script to obtain the conversion factor between physical fluxes and luminosities, must be run within sherpa.

- code/generate_latex_table_code.py: a script to make the latex code for Table 1, 2, and 3 in the paper.

data/: includes all the information of individual LRDs and outputs of CSTACK runs.

- data/all_lrds_final.fits: a file containing all the information we used and derived for our LRD sample. More detailed description of each column can be found in data/explanation.txt.

- data/cstack_output: contains the raw output from CSTACK

- data/distriubtions: contains the derived posterior distribution of all the LRDs and the stack.


More information coming. 

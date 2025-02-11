# nptool_ASGARD
Heads-up: nptool is needed to run this simulation. Please see https://nptool.in2p3.fr/ for details on installation and basic information.

This contains nptool files dedicated for ASGARD (and STARK Jr) simulations. However, you will also find IDATEN-related files as well.
The entire NPTool software is included in this repository, with additional library, detector, input files for ASGARD and STARK Jr detectors. 
For those with NPTool already installed, only the new files will be useful.

A sample detector setup can be found in Inputs/DetectorConfiguration/ASGARD_STARKjr.detector or ASGARD_full_STARKjr.detector.
A sample reaction file (for Coulomb excitation) is found in Inputs/EventGenerator/40Ar_coulex_197Au.reaction.

Sample commands to run npsimulation with ASGARD + STARK Jr. configuration are found in NPSimulation/asgard_starkjr_commands.txt.
A sample command is:
$NPSimulation -D ASGARD_STARKjr.detector -E 40Ar_coulex_197Au_10MeVu_pencil.reaction -B 152Eu_source.mac -O asgard_starkjr_40Ar_10MeVu_pencil

For the reactions of interest, SRIM output files for the beam/target nuclei are needed in .txt format. They should be placed in the Outputs/Analysis directory.

Analysis of output files is done with Outputs/Analysis/asgard_starkjr_analysis.C or asgard_starkjr_10MeVu_analysis.C
It will produce several plots: 
1. raw detected particle energies vs z-position on STARK Jr., corrected energies based on ideal kinematics and calculating back energy lost through the target.
2. 2D energy matrices of Doppler-corrected gamma-ray energies vs cos(theta_pg), a) for beam and b) for target nuclei
3. 1D energy spectra of Doppler-corrected gamma-ray energies.
 
One needs to specify the NPTool output file path in the script and run:
$root asgard_starkjr_10MeVu_analysis.C

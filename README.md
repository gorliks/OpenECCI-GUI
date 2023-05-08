# OpenECCI
GUI-based software for controlled ECCI (Electron Channelling Contrast Imaging) analysis of crystal defects in a scanning electron microscope (SEM).

The analysis relies on electron channeling pattern (ECP) from a reference sample, typically single-crystal Si.
This reference imaging is used to calibrate the SEM stage rotation, in particular the sample stage tilts (tilt and rotation) to match the simulated pattern. Optionally, it is possible to check the electron backscatter diffraction (EBSD)
pattern from the reference sample (typically single-crystal Si as well).

Next, by loading the master pattern generated from EMsoft for the characterised sample
(in the example, a polycristalline Austenitic steel sample with face centre cubic lattice), 
and its EBSD map. For each grain in the sample, a corresponding RKP (Reflected Kikuchi Pattern) can be simulated based on the EBSD Euler angles and the calibration parameters. Stage tilt and rotations can now be calculated in order to bring the selected
Kikuchi band edge or "two beam" conditions to the RKP pattern centre.

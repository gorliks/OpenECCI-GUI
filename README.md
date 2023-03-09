# OpenECCI
GUI-based software for dislocations imaging 
in a scanning electron microscope (SEM).

The analysis relies on electron channeling contrast imaging (ECCI) of a 
reference sample, typically single-crystal Si.
This reference imaging is used to calibrate the SEM, in particular the
sample stage tilts (tilt and rotation) to match the simulated pattern.
Optionally, it is possible to check the electron backscatter diffraction (EBSD)
pattern from the reference sample (typically single-crystal Si as well).

Next, by loading the master pattern for the characterised sample
(in the example, polycristalline Fe face-cubic-symmetric crystal), 
and its EBSD analysis, for each grain in the sample, 
a corresponding ECP (electron channelling pattern) is calculated, 
using the calibration parameters.
Stage tilt and rotations can now be calculated in order to bring the selected
Kikuchi band in the centre of the pattern, which corresponds to a two-vector
condition, when the product of the dislocation Burger vector and the selected 
lattice vector is zero.
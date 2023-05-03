# OpenECCI
GUI-based software for dislocations imaging 
in a scanning electron microscope (SEM).

The analysis relies on electron channeling pattern (ECP) of a 
reference sample, typically single-crystal Si.
This reference imaging is used to calibrate the SEM, in particular the
sample stage tilts (tilt and rotation) to match the simulated pattern.
Optionally, it is possible to check the electron backscatter diffraction (EBSD)
pattern from the reference sample (typically single-crystal Si as well).

Next, by loading the master pattern generated from EMsoft for the characterised sample
(in the example, a polycristalline Fe Austenite face centre cubic crystal), 
and its EBSD analysis, for each grain in the sample, 
a corresponding RKP (Reflected Kikuchi Pattern) is simulated, 
using the calibration parameters.
Stage tilt and rotations can now be calculated in order to bring the selected
Kikuchi band edge in the centre of the pattern, which corresponds to a "two beam"
condition. 

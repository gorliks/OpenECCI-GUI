import numpy as np
import matplotlib.pyplot as plt
import os, glob
import sys
sys.path.insert(0, '../src/')
import utils

#from .OpenECCI.src import utils


import importlib
importlib.reload(utils)

import electron_diffraction


path_to_data = r'C:\Users\sergeyg\Github\OpenECCI\data'
master_pattern_file = r'Si-master-20kv.h5'
ctf_file = r'20kv_26nA_15mm WD_4x4 binning Si Map Data 2.ctf'

measured_ref_ECP_file_name = os.path.join(path_to_data, 'Si_ECP_001.tif')
measured_ref_EBSD_file_name = os.path.join(path_to_data, 'Si_pattern.tiff')

file_name_MP = os.path.join(path_to_data, master_pattern_file)
print(file_name_MP)
file_name_ctf = os.path.join(path_to_data, ctf_file)
print(file_name_ctf)

########################################################################################################################################
ecp_reference = electron_diffraction.Kikuchi(mode='ECP')
ecp_reference.update_settings(energy=20,
                              pc_x=0.5,
                              pc_y=0.5,
                              pc_z=4.0,
                              pixel_size=10,
                              binning=1,
                              detector_tilt=-90,
                              sample_tilt=0,
                              projection='lambert',
                              hemispheres='both',
                              detector_shape=(884, 1024),
                              convention=None)

ecp_reference.load_master_pattern(path_to_master_pattern=file_name_MP)
ecp_reference.load_xmap(file_name=file_name_ctf)
print(f'ECP point Euler angle = {ecp_reference.point_Eulers}')

measured_ref_ECP = utils.load_image(measured_ref_ECP_file_name)
measured_ref_ECP = measured_ref_ECP[:884, :]
simulated_ECP = ecp_reference.calculate_diffraction_pattern(tilt_x=0,
                                                            tilt_y=0)
difference_ECP = utils.calculate_difference(image1=measured_ref_ECP,
                                            image2=simulated_ECP)

########################################################################################################################################
ebsd_reference = electron_diffraction.Kikuchi(mode='EBSD')
ebsd_reference.update_settings(energy=20,
                               pc_x=0.4895443444173796,
                               pc_y=0.62040462473423286,
                               pc_z=0.49422998235067966,
                               pixel_size=59.2,
                               binning=4,
                               detector_tilt=0.0,
                               sample_tilt=70.0,
                               projection='lambert',
                               hemispheres='both',
                               detector_shape=(120, 160),
                               convention='oxford')
ebsd_reference.load_master_pattern(path_to_master_pattern=file_name_MP)
ebsd_reference.load_xmap(file_name=file_name_ctf)
print(f'EBSD point Euler angle = {ebsd_reference.point_Eulers}')

measured_ref_EBSD = utils.load_image(measured_ref_EBSD_file_name)
measured_ref_EBSD = np.flipud(measured_ref_EBSD)
simulated_EBSD = ebsd_reference.calculate_diffraction_pattern(tilt_x=0,
                                                              tilt_y=0)
difference_EBSD = utils.calculate_difference(image1=measured_ref_EBSD,
                                             image2=simulated_EBSD)

########################################################################################################################################

fig = plt.figure(figsize=(10, 5))
ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)
ax6 = fig.add_subplot(236)

ax1.imshow(measured_ref_ECP, cmap='gray')
ax2.imshow(simulated_ECP, cmap='gray')
ax3.imshow(difference_ECP, cmap='bwr')
###
ax4.imshow(measured_ref_EBSD, cmap='gray')
ax5.imshow(simulated_EBSD, cmap='gray')
ax6.imshow(difference_EBSD, cmap='bwr')

plt.show()

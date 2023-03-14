import qtdesigner_files.main_gui as gui_main
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtWidgets import QFileDialog

# import qimage2ndarray
# from importlib import reload  # Python 3.4+
# from dataclasses import dataclass

import sys, time, os, glob
import numpy as np
import copy
import importlib

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as _FigureCanvas
from matplotlib.backends.backend_qt5agg import (
    NavigationToolbar2QT as _NavigationToolbar)

import electron_diffraction
importlib.reload(electron_diffraction)

import utils
importlib.reload(utils)

from diffpy.structure import Atom, Lattice, Structure
from orix.crystal_map import CrystalMap, Phase, PhaseList
from orix.quaternion import Orientation, Rotation, symmetry
from orix.vector import Vector3d, AxAngle, Miller





class GUIMainWindow(gui_main.Ui_MainWindow, QtWidgets.QMainWindow):
    def __init__(self):
        super(GUIMainWindow, self).__init__()
        self.setupUi(self)
        self.path_to_ref_master_pattern = None
        self.path_to_ref_ctf_file = None
        self.path_to_sample_master_pattern = None
        self.path_to_sample_ctf_file = None

        self.measured_ref_ECP = None
        self.measured_ref_ECP_stored = None
        self.simulated_ECP = None
        self.simulated_ECP_stored = None
        self.difference = None

        self.measured_ref_EBSD = None
        self.measured_ref_EBSD_stored = None
        self.simulated_EBSD = None
        self.difference_EBSD = None

        self._abort_clicked_status = False
        self.tilt_x = 0.0
        self.tilt_y = 0.0
        self.calibrated_tilt_x = 0
        self.calibrated_tilt_y = 0
        self.mac = 0.0
        self.mac_max = 0.0

        self.ecp_reference  = electron_diffraction.Kikuchi(mode='ECP')
        self.ebsd_reference = electron_diffraction.Kikuchi(mode='EBSD')
        ###############################################################
        self.ebsd_sample    = electron_diffraction.Kikuchi(mode='EBSD')
        self.ecp_sample     = electron_diffraction.Kikuchi(mode='ECP')
        self.Euler1 = 0
        self.Euler2 = 0
        self.Euler3 = 0
        self.wide_angle_EBSD_pattern = None
        self.wide_angle_ECP_pattern = None
        self.lines = None
        self.zone_axes = None
        self.zone_axes_labels = None

        self.horizontalScrollBar_tilt_y.setValue(50)
        self.horizontalScrollBar_tilt_x.setValue(50)
        self.doubleSpinBox_tilt_y.setValue(0)
        self.doubleSpinBox_tilt_x.setValue(0)
        self.comboBox_angle_step.setCurrentText('0.1')
        self.comboBox_angle_step_automatic.setCurrentText('0.1')
        number_of_steps = 2*self.doubleSpinBox_angle_range.value() / float(self.comboBox_angle_step.currentText())
        self.horizontalScrollBar_tilt_x.setMaximum(int(number_of_steps))
        self.horizontalScrollBar_tilt_y.setMaximum(int(number_of_steps))

        self.setStyleSheet("""QLCDNumber {
        border: 1px solid lightgray;
        border-radius: 5px;
        background-color: rgb(0, 0, 0);
        color: rgb(255, 0, 0)
        }""")
        self.setStyleSheet("""QPushButton {
        border: 1px solid lightgray;
        border-radius: 5px;
        background-color: #e3e3e3;
        }""")

        self.setup_connections()
        self.initialise_image_frames()


    def setup_connections(self):
        self.label_messages.setText('Starting up...')
        self.doubleSpinBox_stage_rotation.setKeyboardTracking(False)
        self.doubleSpinBox_stage_tilt.setKeyboardTracking(False)
        #
        self.pushButton_abort_automatic_calibration.clicked.connect(lambda: self._abort_clicked())
        #
        self.horizontalScrollBar_tilt_y.valueChanged.connect(lambda: self._set_tilt(selector=2, plot=True))
        self.horizontalScrollBar_tilt_x.valueChanged.connect(lambda: self._set_tilt(selector=1, plot=True))
        # self.horizontalScrollBar_tilt_y.sliderReleased.connect(lambda: self._set_tilt3(selector=2, plot=True))
        # self.horizontalScrollBar_tilt_x.sliderReleased.connect(lambda: self._set_tilt3(selector=1, plot=True))
        #
        self.doubleSpinBox_tilt_y.editingFinished.connect(lambda: self._set_tilt2(selector=2, plot=True))
        self.doubleSpinBox_tilt_x.editingFinished.connect(lambda: self._set_tilt2(selector=1, plot=True))
        self.comboBox_angle_step.currentTextChanged.connect(lambda: self._change_angle_step())
        self.doubleSpinBox_angle_range.editingFinished.connect(lambda: self._change_angle_step())
        #
        self.pushButton_open_file_ref_ECCI_measurement.clicked.connect(lambda: self._open_ref_ECCI_measurement_file())
        self.pushButton_crop_meas_ECP_X.clicked.connect(lambda: self._crop_measured_ref_ECP(mode='X'))
        self.pushButton_crop_meas_ECP_Y.clicked.connect(lambda: self._crop_measured_ref_ECP(mode='Y'))
        self.pushButton_crop_meas_ECP_restore.clicked.connect(lambda: self._restore_loaded_pattern())
        self.pushButton_load_ECP_master_pattern.clicked.connect(lambda: self._load_ref_master_pattern())
        self.pushButton_set_ref_ECP_detector.clicked.connect(lambda: self._update_ecp_ref_settings())
        self.pushButton_load_ref_ctf_file.clicked.connect(lambda: self._load_reference_ctf_file())
        self.pushButton_ref_ECP_display.clicked.connect(lambda: self.calculate_simulated_ECP_pattern())
        self.pushButton_run_automatic_calibration.clicked.connect(lambda: self.run_automatic_calibration())
        #
        self.pushButton_open_file_ref_EBSD_measurement.clicked.connect(lambda: self._open_ref_EBSD_measurement_file())
        self.pushButton_set_ref_EBSD_detector.clicked.connect(lambda: self._update_ebsd_ref_settings())
        self.pushButton_calculate_EBSD_pattern_and_display.clicked.connect(lambda:
                                                                           self.calculate_simulated_EBSD_pattern())
        self.pushButton_open_stack_ref_EBSD_measurement.clicked.connect(lambda:
                                                                        self._open_ref_EBSD_measurement_stack())
        self.pushButton_ref_EBSD_generate_dictionary_index.clicked.connect(lambda: self.generate_reference_DI())
        self.pushButton_display_selected_EBSD_pattern.clicked.connect(lambda: self.display_selected_EBSD_pattern())
        self.pushButton_display_calculated_EBSD_for_Eulers.clicked.connect(lambda: self.display_simulated_EBSD_for_Eulers())
        #
        self.pushButton_load_sample_master_pattern.clicked.connect(lambda: self._load_sample_master_pattern())
        self.pushButton_load_sample_ctf_file.clicked.connect(lambda: self._load_sample_ctf_file())
        #
        self.pushButton_set_sample_ECP_detector.clicked.connect(lambda: self._update_ecp_sample_settings())
        self.pushButton_set_sample_EBSD_detector.clicked.connect(lambda: self._update_ebsd_sample_settings())
        self.pushButton_display_for_manual_Eulers.clicked.connect(lambda: self.display_EBSD_ECP_for_Eulers())
        self.checkBox_plot_indexed.clicked.connect(lambda: self.plot_wide_angle_EBSD(Euler_angles=[self.Euler1,
                                                                                                   self.Euler2,
                                                                                                   self.Euler3]))
        self.doubleSpinBox_stage_rotation.valueChanged.connect(lambda: self.plot_wide_angle_EBSD(Euler_angles=[self.Euler1,
                                                                                                               self.Euler2,
                                                                                                               self.Euler3]))
        self.doubleSpinBox_stage_tilt.valueChanged.connect(lambda: self.plot_wide_angle_EBSD(Euler_angles=[self.Euler1,
                                                                                                           self.Euler2,
                                                                                                           self.Euler3]))
        self.doubleSpinBox_stage_rotation.valueChanged.connect(lambda: self.plot_wide_angle_ECP(Euler_angles=[self.Euler1,
                                                                                                              self.Euler2,
                                                                                                              self.Euler3]))
        self.doubleSpinBox_stage_tilt.valueChanged.connect(lambda: self.plot_wide_angle_ECP(Euler_angles=[self.Euler1,
                                                                                                          self.Euler2,
                                                                                                          self.Euler3]))



    def initialise_image_frames(self):
        self.figure_ECCI_exp = plt.figure(10)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_ECCI_exp = _FigureCanvas(self.figure_ECCI_exp)
        self.toolbar_ECCI_exp = _NavigationToolbar(self.canvas_ECCI_exp, self)
        #
        self.label_image_ref_ECCI_measurement.setLayout(QtWidgets.QVBoxLayout())
        self.label_image_ref_ECCI_measurement.layout().addWidget(self.toolbar_ECCI_exp)
        self.label_image_ref_ECCI_measurement.layout().addWidget(self.canvas_ECCI_exp)

        self.figures = {1 : {'fig' : self.figure_ECCI_exp, 'canvas': self.canvas_ECCI_exp, 'toolbar': self.toolbar_ECCI_exp}  }
        ################################################################################################

        self.figure_ECCI_sim = plt.figure(11)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_ECCI_sim = _FigureCanvas(self.figure_ECCI_sim)
        self.toolbar_ECCI_sim = _NavigationToolbar(self.canvas_ECCI_sim, self)
        #
        self.label_image_ref_ECCI_simulation.setLayout(QtWidgets.QVBoxLayout())
        self.label_image_ref_ECCI_simulation.layout().addWidget(self.toolbar_ECCI_sim)
        self.label_image_ref_ECCI_simulation.layout().addWidget(self.canvas_ECCI_sim)

        self.figures[2] = {'fig' : self.figure_ECCI_sim, 'canvas': self.canvas_ECCI_sim, 'toolbar': self.toolbar_ECCI_sim}
        ################################################################################################

        self.figure_ECCI_diff = plt.figure(12)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_ECCI_diff = _FigureCanvas(self.figure_ECCI_diff)
        self.toolbar_ECCI_diff = _NavigationToolbar(self.canvas_ECCI_diff, self)
        #
        self.label_image_ref_ECCI_difference.setLayout(QtWidgets.QVBoxLayout())
        self.label_image_ref_ECCI_difference.layout().addWidget(self.toolbar_ECCI_diff)
        self.label_image_ref_ECCI_difference.layout().addWidget(self.canvas_ECCI_diff)

        self.figures[3] = {'fig' : self.figure_ECCI_diff, 'canvas': self.canvas_ECCI_diff, 'toolbar': self.toolbar_ECCI_diff}
        ################################################################################################

        self.figure_EBSD_exp = plt.figure(13)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_EBSD_exp = _FigureCanvas(self.figure_EBSD_exp)
        self.toolbar_EBSD_exp = _NavigationToolbar(self.canvas_EBSD_exp, self)
        #
        self.label_image_ref_EBSD_measurement.setLayout(QtWidgets.QVBoxLayout())
        self.label_image_ref_EBSD_measurement.layout().addWidget(self.toolbar_EBSD_exp)
        self.label_image_ref_EBSD_measurement.layout().addWidget(self.canvas_EBSD_exp)

        self.figures[4] = {'fig' : self.figure_EBSD_exp, 'canvas': self.canvas_EBSD_exp, 'toolbar': self.toolbar_EBSD_exp}
        ################################################################################################

        self.figure_EBSD_sim = plt.figure(14)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_EBSD_sim = _FigureCanvas(self.figure_EBSD_sim)
        self.toolbar_EBSD_sim = _NavigationToolbar(self.canvas_EBSD_sim, self)
        #
        self.label_image_ref_EBSD_simulation.setLayout(QtWidgets.QVBoxLayout())
        self.label_image_ref_EBSD_simulation.layout().addWidget(self.toolbar_EBSD_sim)
        self.label_image_ref_EBSD_simulation.layout().addWidget(self.canvas_EBSD_sim)

        self.figures[5] = {'fig' : self.figure_EBSD_sim, 'canvas': self.canvas_EBSD_sim, 'toolbar': self.toolbar_EBSD_sim}
        ################################################################################################

        self.figure_EBSD_diff = plt.figure(15)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_EBSD_diff = _FigureCanvas(self.figure_EBSD_diff)
        self.toolbar_EBSD_diff = _NavigationToolbar(self.canvas_EBSD_diff, self)
        #
        self.label_image_ref_EBSD_difference.setLayout(QtWidgets.QVBoxLayout())
        self.label_image_ref_EBSD_difference.layout().addWidget(self.toolbar_EBSD_diff)
        self.label_image_ref_EBSD_difference.layout().addWidget(self.canvas_EBSD_diff)

        self.figures[6] = {'fig' : self.figure_EBSD_diff, 'canvas': self.canvas_EBSD_diff, 'toolbar': self.toolbar_EBSD_diff}
        ################################################################################################

        self.figure_EBSD_stack = plt.figure(16)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_EBSD_stack = _FigureCanvas(self.figure_EBSD_stack)
        self.toolbar_EBSD_stack = _NavigationToolbar(self.canvas_EBSD_stack, self)
        #
        self.label_image_ref_EBSD_stack_measurement.setLayout(QtWidgets.QVBoxLayout())
        self.label_image_ref_EBSD_stack_measurement.layout().addWidget(self.toolbar_EBSD_stack)
        self.label_image_ref_EBSD_stack_measurement.layout().addWidget(self.canvas_EBSD_stack)

        self.figures[7] = {'fig' : self.figure_EBSD_stack, 'canvas': self.canvas_EBSD_stack, 'toolbar': self.toolbar_EBSD_stack}
        ################################################################################################

        self.figure_EBSD_Euler = plt.figure(17)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_EBSD_Euler = _FigureCanvas(self.figure_EBSD_Euler)
        self.toolbar_EBSD_Euler = _NavigationToolbar(self.canvas_EBSD_Euler, self)
        #
        self.label_image_ref_EBSD_simulation_Euler_angles.setLayout(QtWidgets.QVBoxLayout())
        self.label_image_ref_EBSD_simulation_Euler_angles.layout().addWidget(self.toolbar_EBSD_Euler)
        self.label_image_ref_EBSD_simulation_Euler_angles.layout().addWidget(self.canvas_EBSD_Euler)

        self.figures[8] = {'fig' : self.figure_EBSD_Euler, 'canvas': self.canvas_EBSD_Euler, 'toolbar': self.toolbar_EBSD_Euler}
        ################################################################################################

        self.figure_EBSD_sample = plt.figure(18)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_EBSD_sample = _FigureCanvas(self.figure_EBSD_sample)
        self.toolbar_EBSD_sample = _NavigationToolbar(self.canvas_EBSD_sample, self)
        #
        self.label_image_EBSD_sample_ctf.setLayout(QtWidgets.QVBoxLayout())
        self.label_image_EBSD_sample_ctf.layout().addWidget(self.toolbar_EBSD_sample)
        self.label_image_EBSD_sample_ctf.layout().addWidget(self.canvas_EBSD_sample)

        self.figures[9] = {'fig' : self.figure_EBSD_sample, 'canvas': self.canvas_EBSD_sample, 'toolbar': self.toolbar_EBSD_sample}
        ################################################################################################

        ################################################################################################
        self.figure_stereo_projection = plt.figure(31)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_stereo_projection = _FigureCanvas(self.figure_stereo_projection)
        self.toolbar_stereo_projection = _NavigationToolbar(self.canvas_stereo_projection, self)
        #
        self.label_image_stereo_projection.setLayout(QtWidgets.QVBoxLayout())
        # self.label_image_stereo_projection.layout().addWidget(self.toolbar_stereo_projection)
        self.label_image_stereo_projection.layout().addWidget(self.canvas_stereo_projection)
        ################################################################################################
        ################################################################################################
        self.figure_EBSD_wide_angle_sim = plt.figure(32)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_EBSD_wide_angle_sim = _FigureCanvas(self.figure_EBSD_wide_angle_sim)
        self.toolbar_EBSD_wide_angle_sim = _NavigationToolbar(self.canvas_EBSD_wide_angle_sim, self)
        #
        self.label_image_sample_EBSD_simulation.setLayout(QtWidgets.QVBoxLayout())
        # self.label_image_sample_EBSD_simulation.layout().addWidget(self.toolbar_EBSD_wide_angle_sim)
        self.label_image_sample_EBSD_simulation.layout().addWidget(self.canvas_EBSD_wide_angle_sim)
        ################################################################################################
        ################################################################################################
        self.figure_ECP_wide_angle_sim = plt.figure(33)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_ECP_wide_angle_sim = _FigureCanvas(self.figure_ECP_wide_angle_sim)
        self.toolbar_ECP_wide_angle_sim = _NavigationToolbar(self.canvas_ECP_wide_angle_sim, self)
        #
        self.label_image_sample_ECP_simulation.setLayout(QtWidgets.QVBoxLayout())
        # self.label_image_sample_ECP_simulation.layout().addWidget(self.toolbar_ECP_wide_angle_sim)
        self.label_image_sample_ECP_simulation.layout().addWidget(self.canvas_ECP_wide_angle_sim)
        ################################################################################################


    def update_display(self,
                       image, mode='ref_ECCI_measurement',
                       centre_marker=True):
        cmap = 'gray'
        if mode=='ref_ECCI_measurement':
            key = 1
        if mode=='ref_ECCI_simulation':
            key = 2
        if mode=='difference':
            key = 3
            cmap = 'bwr'
        if mode=='ref_EBSD_measurement':
            key = 4
        if mode=='ref_EBSD_simulation':
            key = 5
        if mode=='EBSD_difference':
            key = 6
            cmap='bwr'
        if mode=='ref_EBSD_stack':
            key = 7
        if mode=='ref_EBSD_simulation_Euler':
            key = 8
        if mode=='sample_EBSD_ctf':
            key = 9

        self.figures[key]['fig'].clear()
        self.figures[key]['fig'].patch.set_facecolor(
                             (240 / 255, 240 / 255, 240 / 255))
        self.ax = self.figures[key]['fig'].add_subplot(111)
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)
        self.ax.imshow(image, cmap=cmap)
        if centre_marker:
            dims = image.shape
            centre_marker_coords = [dims[0]/2, dims[1]/2]
            self.ax.plot(centre_marker_coords[1],
                         centre_marker_coords[0],
                         '+', c='red', markersize=15, markeredgewidth=3)
        self.figures[key]['canvas'].draw()

        if mode=='sample_EBSD_ctf':
            def on_click(event):
                coords = []
                coords.append(event.ydata)
                coords.append(event.xdata)
                try:
                    coords = np.flip(coords[-2:], axis=0)
                    x_clicked = coords[0]
                    y_clicked = coords[1]
                except:
                    x_clicked = 0
                    y_clicked = 0

                [Eu1, Eu2, Eu3] = np.rad2deg(
                    Rotation.to_euler(
                        self.ebsd_sample.xmap[int(y_clicked),
                                              int(x_clicked)].orientations))[0]

                print(' - - - - - > CLICKED:', x_clicked, y_clicked, Eu1, Eu2, Eu3)
                self.lcdNumber_Euler1.display(Eu1)
                self.lcdNumber_Euler2.display(Eu2)
                self.lcdNumber_Euler3.display(Eu3)
                self.Euler1 = Eu1
                self.Euler2 = Eu2
                self.Euler3 = Eu3
                self.doubleSpinBox_Euler_1.setValue(Eu1)
                self.doubleSpinBox_Euler_2.setValue(Eu2)
                self.doubleSpinBox_Euler_3.setValue(Eu3)

                self.plot_stereo_projection(Euler_angles=[Eu1, Eu2, Eu3])
                self.plot_wide_angle_EBSD(Euler_angles=[Eu1, Eu2, Eu3])
                self.plot_wide_angle_ECP(Euler_angles=[Eu1, Eu2, Eu3])

            self.figures[key]['fig'].canvas.mpl_connect("button_press_event",
                                                        on_click)
            self.blitted_cursor = utils.BlittedCursor(self.ax)
            self.figures[key]['fig'].canvas.mpl_connect('motion_notify_event',
                                                        self.blitted_cursor.on_mouse_move)




    def plot_wide_angle_EBSD(self, Euler_angles=[0, 0, 0]):
        self._update_ebsd_sample_settings()
        Euler_angles = np.radians(Euler_angles)
        self.wide_angle_EBSD_pattern = \
            self.ebsd_sample.calculate_diffraction_pattern(tilt_x=self.doubleSpinBox_tilt_x_calibrated.value(),
                                                           tilt_y=self.doubleSpinBox_tilt_y_calibrated.value(),
                                                           stage_rotation=self.doubleSpinBox_stage_rotation.value(),
                                                           stage_tilt=self.doubleSpinBox_stage_tilt.value(),
                                                           Eulers=Euler_angles)
        self.figure_EBSD_wide_angle_sim.clear()
        self.figure_EBSD_wide_angle_sim.patch.set_facecolor(
                                 (240 / 255, 240 / 255, 240 / 255))
        self.ax = self.figure_EBSD_wide_angle_sim.add_subplot(111)
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)
        self.ax.imshow(self.wide_angle_EBSD_pattern, cmap='gray')
        if self.checkBox_plot_indexed.isChecked():
            self.lines, self.zone_axes, self.zone_axes_labels = \
                self.ebsd_sample.get_indexed_kikuchi(tilt_x=self.doubleSpinBox_tilt_x_calibrated.value(),
                                                     tilt_y=self.doubleSpinBox_tilt_y_calibrated.value(),
                                                     Euler_angles=Euler_angles,
                                                     stage_rotation=self.doubleSpinBox_stage_rotation.value(),
                                                     stage_tilt=self.doubleSpinBox_stage_tilt.value())
            self.ax.add_collection(self.lines)
            self.ax.add_collection(self.zone_axes)
            for label in self.zone_axes_labels:
                self.ax.add_artist(label)

        self.canvas_EBSD_wide_angle_sim.draw()



    def plot_wide_angle_ECP(self, Euler_angles=[0, 0, 0]):
        self._update_ecp_sample_settings()
        Euler_angles = np.radians(Euler_angles)
        self.wide_angle_ECP_pattern = \
            self.ecp_sample.calculate_diffraction_pattern(tilt_x=self.doubleSpinBox_tilt_x_calibrated.value(),
                                                          tilt_y=self.doubleSpinBox_tilt_y_calibrated.value(),
                                                          stage_rotation=self.doubleSpinBox_stage_rotation.value(),
                                                          stage_tilt=self.doubleSpinBox_stage_tilt.value(),
                                                          Eulers=Euler_angles)
        self.figure_ECP_wide_angle_sim.clear()
        self.figure_ECP_wide_angle_sim.patch.set_facecolor(
                                (240 / 255, 240 / 255, 240 / 255))
        self.ax = self.figure_ECP_wide_angle_sim.add_subplot(111)
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)
        self.ax.imshow(self.wide_angle_ECP_pattern, cmap='gray')

        dims = self.wide_angle_ECP_pattern.shape
        centre_marker_coords = [dims[0] / 2, dims[1] / 2]
        self.ax.plot(centre_marker_coords[1],
                     centre_marker_coords[0],
                     '+', c='red', markersize=15, markeredgewidth=3)

        self.canvas_ECP_wide_angle_sim.draw()

        def on_click(event):
            coords = []
            coords.append(event.ydata)
            coords.append(event.xdata)
            try:
                x_pos = coords[-2]
                y_pos = coords[-1]
            except:
                x_pos = 0
                y_pos = 0

            [PCx_bkp, PCy_bkp, PCz_bkp] = self.ecp_sample.detector.pc[0]
            [Ny, Nx] = self.ecp_sample.detector.shape
            px_size_bkp =  self.ecp_sample.detector.px_size
            binning_bkp =  self.ecp_sample.detector.binning

            # Calculated the phycial distances on the detector from selected pixel to the projection centre, calculated in micrometer/um
            [coord_x, coord_y] = [x_pos, y_pos]
            distance_x = (coord_x - PCx_bkp * Nx) * px_size_bkp * binning_bkp
            distance_y = (PCy_bkp * Ny - coord_y) * px_size_bkp * binning_bkp
            distance_l = PCz_bkp * Ny * px_size_bkp * binning_bkp

            # Calcuated the azimuthal and polar angle of the selected pixel in the detector frame
            if distance_x == 0 and distance_y == 0:
                azi_bkp = 0
            elif distance_x > 0 and distance_y == 0:
                azi_bkp = 0
            elif distance_x == 0 and distance_y > 0:
                azi_bkp = np.pi / 2
            elif distance_x < 0 and distance_y == 0:
                azi_bkp = np.pi
            elif distance_x == 0 and distance_y < 0:
                azi_bkp = np.pi * 3 / 2
            elif distance_x < 0:
                azi_bkp = np.arctan(distance_y / distance_x) + np.pi
            elif distance_x > 0 and distance_y < 0:
                azi_bkp = np.arctan(distance_y / distance_x) + np.pi * 2
            else:
                azi_bkp = np.arctan(distance_y / distance_x)

            polar_bkp = np.arctan(np.sqrt(distance_x ** 2 + distance_y ** 2) / distance_l)

            # self.doubleSpinBox_stage_rotation.setValue(np.degrees(azi_bkp)-90)
            # self.doubleSpinBox_stage_tilt.setValue(np.degrees(polar_bkp))
            self.label_messages.setText(f"Pixel position {int(x_pos)}, {int(y_pos)}\n "
                                        f"Physical distance {round(distance_x,2), round(distance_y,2), round(distance_l,2)}um \n"
                                        f"Stage Rot {round(np.degrees(azi_bkp)-90,2)}\N{DEGREE SIGN}, "
                                        f"Stage tilt {round(np.degrees(polar_bkp),2)}\N{DEGREE SIGN}")

        self.figure_ECP_wide_angle_sim.canvas.mpl_connect("button_press_event",
                                                          on_click)
        self.blitted_cursor2 = utils.BlittedCursor(self.ax)
        self.figure_ECP_wide_angle_sim.canvas.mpl_connect('motion_notify_event',
                                                          self.blitted_cursor2.on_mouse_move)



    def plot_stereo_projection(self, Euler_angles = [0,0,0]):
        print('stereo projection: ', Euler_angles)

        n = int(90 / 10)  # Degree / net resolution
        steps = 500
        kwargs = dict(linewidth=0.25, color="green")
        polar = np.linspace(0, 0.5 * np.pi, num=n)
        v_right = Vector3d.from_polar(azimuth=np.zeros(n), polar=polar)
        v_left = Vector3d.from_polar(azimuth=np.ones(n) * np.pi, polar=polar)
        v010 = Vector3d.zero(shape=(n,))
        v010.y = 1
        v010_opposite = -v010
        #################################################################
        self.figure_stereo_projection.clear()
        self.ax = self.figure_stereo_projection.subplots(subplot_kw=dict(projection="stereographic"))
        self.ax.draw_circle(v_right, steps=steps, **kwargs)
        self.ax.draw_circle(v_left,  steps=steps, **kwargs)
        self.ax.draw_circle(v010, opening_angle=polar, steps=steps, **kwargs)
        self.ax.draw_circle(v010_opposite, opening_angle=polar, steps=steps, **kwargs)
        for label, azimuth in zip(["", "", "", ""], np.array([0, 0.5, 1, 1.5]) * np.pi):
            self.ax.text(azimuth, 0.5 * np.pi, s=label, c="C1")

        cubic = Phase(point_group="m-3m", structure=Structure())
        #m8 = Miller(uvw=[[1, 0, 0]], phase=cubic)
        millers = self.plainTextEdit_m8.toPlainText()
        try:
            millers = np.fromstring(millers, dtype=int, sep=',')
            if len(millers)>3: #more than 1 vector
                millers = np.reshape(millers, (3,-1))
        except:
            millers = [[1, 0, 0], [1,1,1]]
        m8 = Miller(uvw=millers, phase=cubic)
        o = Orientation.from_euler(np.deg2rad(Euler_angles))
        m9, idx = m8.symmetrise(unique=True, return_index=True)
        m9 = ~o * m9
        m9[idx == 0].scatter(hemisphere="upper", c=f"C0",
                             figure=self.figure_stereo_projection,
                             axes_labels=["X", "Y", "Z"], )
        # renderer = self.figure_stereo_projection.canvas.renderer
        self.canvas_stereo_projection.draw()


    def display_EBSD_ECP_for_Eulers(self):
        Eu1 = self.doubleSpinBox_Euler_1.value()
        Eu2 = self.doubleSpinBox_Euler_2.value()
        Eu3 = self.doubleSpinBox_Euler_3.value()
        Eulers = np.array([Eu1, Eu2, Eu3])
        self.plot_wide_angle_EBSD(Euler_angles=Eulers)
        self.plot_wide_angle_ECP(Euler_angles=Eulers)



    ####################################################################################################################

    def _update_ecp_ref_settings(self):
        if self.measured_ref_ECP is not None:
            energy = self.doubleSpinBox_ref_ECP_energy.value()
            projection = self.comboBox_ref_ECP_projection.currentText()
            hemispheres = self.comboBox_ref_ECP_hemispheres.currentText()
            shape = self.measured_ref_ECP.shape
            pc_x = self.doubleSpinBox_ref_ECP_pc_x.value()
            pc_y = self.doubleSpinBox_ref_ECP_pc_y.value()
            pc_z = self.doubleSpinBox_ref_ECP_pc_z.value()
            pixel_size = self.doubleSpinBox_ref_ECP_pixel_size.value()
            detector_tilt = self.doubleSpinBox_ref_ECP_detector_tilt.value()
            sample_tilt = self.doubleSpinBox_ref_ECP_sample_tilt.value()
            binning = self.spinBox_ref_ECP_binning.value()
            self.ecp_reference.update_settings(energy=energy,
                                               pc_x = pc_x,
                                               pc_y = pc_y,
                                               pc_z = pc_z,
                                               pixel_size = pixel_size,
                                               binning = binning,
                                               detector_tilt = detector_tilt,
                                               sample_tilt = sample_tilt,
                                               projection = projection,
                                               hemispheres = hemispheres,
                                               detector_shape = shape)
            self.label_messages.setText('ref ECP detector settings uploaded')
            print(energy, pc_x, pc_y, pc_z, pixel_size, binning, detector_tilt, sample_tilt, projection, hemispheres, shape)
        else:
            self.label_messages.setText('First, load a measured ECP/ECCI pattern to define the number of pixels in the pattern')


    def _update_ebsd_ref_settings(self):
        if self.measured_ref_EBSD is not None:
            energy = self.doubleSpinBox_ref_EBSD_energy.value()
            projection = self.comboBox_ref_EBSD_projection.currentText()
            hemispheres = self.comboBox_ref_EBSD_hemispheres.currentText()
            shape = self.measured_ref_EBSD.shape
            pc_x = self.doubleSpinBox_ref_EBSD_pc_x.value()
            pc_y = self.doubleSpinBox_ref_EBSD_pc_y.value()
            pc_z = self.doubleSpinBox_ref_EBSD_pc_z.value()
            pixel_size = self.doubleSpinBox_ref_EBSD_pixel_size.value()
            detector_tilt = self.doubleSpinBox_ref_EBSD_detector_tilt.value()
            sample_tilt = self.doubleSpinBox_ref_EBSD_sample_tilt.value()
            binning = self.spinBox_ref_EBSD_binning.value()

            convention = self.comboBox_ref_EBSD_convention.currentText()

            self.ebsd_reference.update_settings(energy=energy,
                                                pc_x = pc_x,
                                                pc_y = pc_y,
                                                pc_z = pc_z,
                                                pixel_size = pixel_size,
                                                binning = binning,
                                                detector_tilt = detector_tilt,
                                                sample_tilt = sample_tilt,
                                                projection = projection,
                                                hemispheres = hemispheres,
                                                detector_shape = shape,
                                                convention=convention)
            self.label_messages.setText('ref EBSD detector settings uploaded')
            if self.path_to_ref_master_pattern is not None:
                self.ebsd_reference.load_master_pattern(path_to_master_pattern=self.path_to_ref_master_pattern)
            if self.path_to_ref_ctf_file is not None:
                self.ebsd_reference.load_xmap(file_name=self.path_to_ref_ctf_file)
            print(energy, pc_x, pc_y, pc_z, pixel_size, binning, detector_tilt, sample_tilt, projection, hemispheres, shape, convention)
        else:
            self.label_messages.setText('First, load a measured EBSD pattern to load metadata and define the number of pixels in the pattern and other settings')

    ####################################################################################################################

    def _update_ecp_sample_settings(self):
        energy = self.doubleSpinBox_sample_energy.value()
        projection = self.comboBox_sample_ECP_projection.currentText()
        hemispheres = self.comboBox_sample_ECP_hemispheres.currentText()
        shape = (self.spinBox_sample_ECP_pattern_pixels_x.value(),
                 self.spinBox_sample_ECP_pattern_pixels_y.value())
        pc_x = self.doubleSpinBox_sample_ECP_pc_x.value()
        pc_y = self.doubleSpinBox_sample_ECP_pc_y.value()
        pc_z = self.doubleSpinBox_sample_ECP_pc_z.value()
        pixel_size = self.doubleSpinBox_sample_ECP_pixel_size.value()
        detector_tilt = self.doubleSpinBox_sample_ECP_detector_tilt.value()
        sample_tilt = self.doubleSpinBox_sample_ECP_sample_tilt.value()
        binning = self.spinBox_sample_ECP_binning.value()
        self.ecp_sample.update_settings(energy=energy,
                                        pc_x = pc_x,
                                        pc_y = pc_y,
                                        pc_z = pc_z,
                                        pixel_size = pixel_size,
                                        binning = binning,
                                        detector_tilt = detector_tilt,
                                        sample_tilt = sample_tilt,
                                        projection = projection,
                                        hemispheres = hemispheres,
                                        detector_shape = shape)
        self.label_messages.setText('sample ECP detector settings uploaded')
        print(energy, pc_x, pc_y, pc_z, pixel_size, binning, detector_tilt, sample_tilt, projection, hemispheres, shape)


    def _update_ebsd_sample_settings(self):
        energy = self.doubleSpinBox_sample_energy.value()
        projection = self.comboBox_sample_EBSD_projection.currentText()
        hemispheres = self.comboBox_sample_EBSD_hemispheres.currentText()
        shape = (self.spinBox_sample_EBSD_pattern_pixels_x.value(),
                 self.spinBox_sample_EBSD_pattern_pixels_y.value())
        pc_x = self.doubleSpinBox_sample_EBSD_pc_x.value()
        pc_y = self.doubleSpinBox_sample_EBSD_pc_y.value()
        pc_z = self.doubleSpinBox_sample_EBSD_pc_z.value()
        pixel_size = self.doubleSpinBox_sample_EBSD_pixel_size.value()
        detector_tilt = self.doubleSpinBox_sample_EBSD_detector_tilt.value()
        sample_tilt = self.doubleSpinBox_sample_EBSD_sample_tilt.value()
        binning = self.spinBox_sample_EBSD_binning.value()

        convention = self.comboBox_sample_EBSD_convention.currentText()

        self.ebsd_sample.update_settings(energy=energy,
                                         pc_x = pc_x,
                                         pc_y = pc_y,
                                         pc_z = pc_z,
                                         pixel_size = pixel_size,
                                         binning = binning,
                                         detector_tilt = detector_tilt,
                                         sample_tilt = sample_tilt,
                                         projection = projection,
                                         hemispheres = hemispheres,
                                         detector_shape = shape,
                                         convention=convention)
        self.label_messages.setText('sample EBSD detector settings uploaded')
        print(energy, pc_x, pc_y, pc_z, pixel_size, binning, detector_tilt, sample_tilt, projection, hemispheres, shape, convention)

    ####################################################################################################################



    def _open_ref_ECCI_measurement_file(self):
        status, file_name, data = self._open_tif_file()
        if status==True:
            self.label_messages.setText('loading tiff file ' + file_name)
            self.measured_ref_ECP = utils.load_image(file_name)
            Ny, Nx = self.measured_ref_ECP.shape
            self.spinBox_ref_ECP_pattern_pixels_x.setValue(Nx)
            self.spinBox_ref_ECP_pattern_pixels_y.setValue(Ny)
            #
            self.update_display(image=self.measured_ref_ECP, mode='ref_ECCI_measurement')


    def _open_ref_EBSD_measurement_file(self):
        status, file_name, data = self._open_tif_file()
        if status==True:
            self.label_messages.setText('loading tiff file ' + file_name)
            self.measured_ref_EBSD = utils.load_image(file_name)

            try:
                metadata = utils.HKL_metadata(file_name)
                self.ebsd_reference.update_settings_from_dict(dict=metadata)

            except Exception as e:
                self.label_messages.setText('Could not read the tif metadata')
                print('Could not read the tiff metadata')

            Ny, Nx = self.measured_ref_EBSD.shape
            self.spinBox_ref_EBSD_pattern_pixels_x.setValue(Nx)
            self.spinBox_ref_EBSD_pattern_pixels_y.setValue(Ny)
            #
            self.update_display(image=self.measured_ref_EBSD, mode='ref_EBSD_measurement')

            self.doubleSpinBox_ref_EBSD_pc_x.setValue(self.ebsd_reference.pc_x)
            self.doubleSpinBox_ref_EBSD_pc_y.setValue(self.ebsd_reference.pc_y)
            self.doubleSpinBox_ref_EBSD_pc_z.setValue(self.ebsd_reference.pc_z)
            self.doubleSpinBox_ref_EBSD_energy.setValue(self.ebsd_reference.energy)
            self.doubleSpinBox_ref_EBSD_detector_tilt.setValue(self.ebsd_reference.detector_tilt)
            self.doubleSpinBox_ref_EBSD_sample_tilt.setValue(self.ebsd_reference.sample_tilt)

            self._update_ebsd_ref_settings()

            if self.path_to_ref_master_pattern is not None:
                self.ebsd_reference.load_master_pattern(path_to_master_pattern=self.path_to_ref_master_pattern)
            if self.path_to_ref_ctf_file is not None:
                self.ebsd_reference.load_xmap(file_name=self.path_to_ref_ctf_file)


    def _open_ref_EBSD_measurement_stack(self):
        """
        Select the file of the (i x j) EBSD measurement of the reference sample
        pattern dimentions are the same as ref_pattern previously loaded (Nx x Ny)
        :return:
        None
        """
        status, file_name = self._open_ebsp_stack()
        if status==True:
            self.label_messages.setText('loading ebsp file ' + file_name)
            i, j, Nx, Ny = \
                self.ebsd_reference.load_stack(file_name)
            self.spinBox_ref_stack_i.setValue(i)
            self.spinBox_ref_stack_j.setValue(j)
            self.spinBox_ref_stack_Nx.setValue(Nx)
            self.spinBox_ref_stack_Ny.setValue(Ny)


    def generate_reference_DI(self):
        """
        Create the dictionary EBSD index for the reference sample
        using the masterpattern previously loaded into the ebsd_reference instance
        :return:
        None
        """
        self.label_messages.setText('This is an extremely memory intensive operation and is not implemented')
        print('This is an extremely memory intensive operation and is not implemented')
        #
        # Uncomment to create the dictionary index
        #self.ebsd_reference.generate_dictionary_index()



    def _load_ref_master_pattern(self):
        status, file_name = self._open_master_pattern()
        if status==True:
            print(file_name)
            self.path_to_ref_master_pattern = file_name

            self.label_messages.setText('ref ECP masterpattern: ' + file_name)
            self.label_ecp_master_pattern_path.setText(file_name)

            self._update_ecp_ref_settings()

            output_ecp = \
                self.ecp_reference.load_master_pattern(path_to_master_pattern=file_name,
                                                       name=self.plainTextEdit_ref_masterpattern_name.toPlainText(),
                                                       space_group=self.spinBox_ref_masterpattern_space_group.value())
            output_ebsd = \
                self.ebsd_reference.load_master_pattern(path_to_master_pattern=file_name,
                                                        name=self.plainTextEdit_ref_masterpattern_name.toPlainText(),
                                                        space_group=self.spinBox_ref_masterpattern_space_group.value())
            self.label_messages.setText(output_ecp + '; ' + output_ebsd)


    def _load_sample_master_pattern(self):
        status, file_name = self._open_master_pattern()
        if status==True:
            print(file_name)
            self.path_to_sample_master_pattern = file_name

            self.label_messages.setText('sample ECP masterpattern: ' + file_name)
            self.label_sample_master_pattern_path.setText(file_name)

            # self._update_ecp_sample_settings()

            output_ecp = \
                self.ecp_sample.load_master_pattern(path_to_master_pattern=file_name,
                                                    name=self.plainTextEdit_sample_masterpattern_name.toPlainText(),
                                                    space_group=self.spinBox_sample_masterpattern_space_group.value())
            output_ebsd = \
                self.ebsd_sample.load_master_pattern(path_to_master_pattern=file_name,
                                                     name=self.plainTextEdit_sample_masterpattern_name.toPlainText(),
                                                     space_group=self.spinBox_sample_masterpattern_space_group.value())
            self.label_messages.setText(output_ecp + '; ' + output_ebsd)


    def _load_reference_ctf_file(self):
        status, file_name = self._open_ctf_file()
        if status==True:
            print(file_name)
            self.path_to_ref_ctf_file = file_name

            self.label_messages.setText('ctf filename: ' + file_name)
            self.label_ref_ctf_file.setText(file_name)

            output_ecp = \
                self.ecp_reference.load_xmap(file_name=file_name)
            output_ebsd = \
                self.ebsd_reference.load_xmap(file_name=file_name)
            self.label_messages.setText(str(output_ecp) + '; ' + str(output_ebsd))


    def _load_sample_ctf_file(self):
        status, file_name = self._open_ctf_file()
        if status==True:
            print(file_name)
            self.label_messages.setText('Loading the ctf file. please hold...')
            self.path_to_sample_ctf_file = file_name

            self.label_messages.setText('ctf filename: ' + file_name)
            self.label_sample_ctf_file.setText(file_name)

            skiprows = self.spinBox_skiprows.value()

            """ - - - -  Create unit cells of the phases - - - - """
            structures = [
                Structure(title="non-index"),
                Structure(title=self.comboBox_sample_structure_type.currentText(),
                          atoms=[Atom(self.plainTextEdit_sample_atom.toPlainText(), [0] * 3)],
                          lattice=Lattice(self.doubleSpinBox_lattice1.value(),
                                          self.doubleSpinBox_lattice2.value(),
                                          self.doubleSpinBox_lattice3.value(),
                                          self.doubleSpinBox_lattice4.value(),
                                          self.doubleSpinBox_lattice5.value(),
                                          self.doubleSpinBox_lattice6.value())  )
            ]
            phase_list = PhaseList(
                names=["non-index", self.comboBox_sample_structure_type.currentText() ],
                point_groups=[None, self.plainTextEdit_sample_point_group.toPlainText()],
                structures=structures,
            )

            output_ecp = \
                self.ecp_sample.load_xmap(file_name=file_name, skiprows=skiprows)
            output_ebsd = \
                self.ebsd_sample.load_xmap_sample(file_name=file_name, skiprows=skiprows,
                                                  correction=+90, phase_list=phase_list)
            self.label_messages.setText(str(output_ecp) + '; ' + str(output_ebsd))

            if self.ebsd_sample.xmap_gb is not None:
                self.update_display(image=self.ebsd_sample.xmap_gb,
                                    mode='sample_EBSD_ctf')
            self.label_messages.setText('ctf file has been loaded.')


    def calculate_simulated_ECP_pattern(self, plot=True):
        try:
            self.simulated_ECP = \
                self.ecp_reference.calculate_diffraction_pattern(tilt_x=self.tilt_x,
                                                                 tilt_y=self.tilt_y)
            self.update_display(image=self.simulated_ECP,
                                mode='ref_ECCI_simulation')

            # if both experimental pattern and simulated pattern exist,
            # calculate their difference and plot if True
            if (self.simulated_ECP is not None) and (self.measured_ref_ECP is not None):
                self.difference = utils.calculate_difference(image1=self.measured_ref_ECP,
                                                             image2=self.simulated_ECP)
                #calculate the similarity of two images, display MAC
                self.mac = utils.modal_assurance_criterion(image1=self.measured_ref_ECP,
                                                           image2=self.simulated_ECP)
                self.doubleSpinBox_mac.setValue(self.mac)

                """  update the best tilt angles for calibration """
                if self.mac > self.mac_max:
                    self.doubleSpinBox_tilt_x_calibrated.setValue(self.tilt_x)
                    self.doubleSpinBox_tilt_y_calibrated.setValue(self.tilt_y)
                    self.calibrated_tilt_x = self.tilt_x
                    self.calibrated_tilt_y = self.tilt_y
                    self.mac_max = self.mac
                    self.doubleSpinBox_mac_max.setValue(self.mac_max)
                    print(f'mac_max = {self.mac_max}')

                if plot==True:
                    self.update_display(image=self.difference,
                                        mode='difference')
        except Exception as e:
            print('could not simulate the pattern', e)


    def calculate_simulated_EBSD_pattern(self, plot=True):
        try:
            self.simulated_EBSD = \
                self.ebsd_reference.calculate_diffraction_pattern(tilt_x=0,
                                                                  tilt_y=0)
            self.update_display(image=self.simulated_EBSD,
                                mode='ref_EBSD_simulation')

            # if both experimental pattern and simulated pattern exist,
            # calculate their difference and plot if True
            if (self.simulated_EBSD is not None) and (self.measured_ref_EBSD is not None):
                self.difference_EBSD = utils.calculate_difference(image1=self.measured_ref_EBSD,
                                                                  image2=self.simulated_EBSD)
                if plot==True:
                    self.update_display(image=self.difference_EBSD,
                                        mode='EBSD_difference')
        except Exception as e:
            print('could not simulate the EBSD pattern', e)


    def display_selected_EBSD_pattern(self):
        if self.ebsd_reference.stack is not None:
            i = self.spinBox_ref_stack_i_to_plot.value()
            j = self.spinBox_ref_stack_j_to_plot.value()
            self.update_display(image=self.ebsd_reference.stack.inav[i,j].data,
                                mode='ref_EBSD_stack')


    def display_simulated_EBSD_for_Eulers(self):
        Euler1 = self.doubleSpinBox_Euler_1.value()
        Euler2 = self.doubleSpinBox_Euler_2.value()
        Euler3 = self.doubleSpinBox_Euler_3.value()
        Eulers = np.array([Euler1, Euler2, Euler3])
        Eulers = np.radians(Eulers)
        simulated = \
            self.ebsd_reference.calculate_diffraction_pattern(tilt_x=0,
                                                              tilt_y=0,
                                                              Eulers=Eulers)
        self.update_display(image=simulated,
                            mode='ref_EBSD_simulation_Euler')



    def run_automatic_calibration(self):
        mac_max = 0
        tilt_x_max = 0
        tilt_y_max = 0
        ########################################################
        tilt_x_0 = self.doubleSpinBox_tilt_x.value()
        tilt_y_0 = self.doubleSpinBox_tilt_y.value()
        half_range = self.doubleSpinBox_angle_range_automatic.value()
        start_x = tilt_x_0 - half_range
        end_x   = tilt_x_0 + half_range
        start_y = tilt_y_0 - half_range
        end_y   = tilt_y_0 + half_range
        step = float(self.comboBox_angle_step_automatic.currentText())
        TILT_X = np.arange(start_x, end_x+step, step)
        TILT_Y = np.arange(start_y, end_y+step, step)
        for tilt_x in TILT_X:
            for tilt_y in TILT_Y:
                self.repaint()
                QtWidgets.QApplication.processEvents()
                if self._abort_clicked_status == True:
                    print('Abort clicked')
                    self._abort_clicked_status = False  # reinitialise back to False
                    return
                self.doubleSpinBox_tilt_x.setValue(tilt_x)
                self.doubleSpinBox_tilt_y.setValue(tilt_y)
                self._set_tilt2(selector=1, plot=False) # update the slider
                self._set_tilt2(selector=2, plot=False) # update the slider
                _plot_ = self.checkBox_plot_while_running.isChecked()
                self.calculate_simulated_ECP_pattern(plot=_plot_)
                if self.mac > mac_max:
                    mac_max = self.mac
                    self.mac_max = self.mac
                    tilt_x_max = tilt_x
                    tilt_y_max = tilt_y
                    self.calibrated_tilt_x = tilt_x_max
                    self.calibrated_tilt_y = tilt_y_max
                    self.doubleSpinBox_tilt_x_calibrated.setValue(tilt_x_max)
                    self.doubleSpinBox_tilt_y_calibrated.setValue(tilt_y_max)
                    self.doubleSpinBox_mac_max.setValue(mac_max)
                self.repaint()
                QtWidgets.QApplication.processEvents()

        # update the calibrated tilt values at the end of the run
        self.calibrated_tilt_x = tilt_x_max
        self.calibrated_tilt_y = tilt_y_max
        self.label_messages.setText('Calibration completed')




    def _crop_measured_ref_ECP(self, mode='X'):
        if self.measured_ref_ECP is not None:
            start = self.spinBox_ref_meas_ECP_crop_start.value()
            end = self.spinBox_ref_meas_ECP_crop_end.value()
            self.measured_ref_ECP_stored = copy.deepcopy(self.measured_ref_ECP)
            if mode=="X":
                self.measured_ref_ECP = self.measured_ref_ECP[:, start : end]
            else:
                self.measured_ref_ECP = self.measured_ref_ECP[start : end, :]
            self.update_display(image=self.measured_ref_ECP, mode='ref_ECCI_measurement')
            Ny, Nx = self.measured_ref_ECP.shape
            self.spinBox_ref_ECP_pattern_pixels_x.setValue(Nx)
            self.spinBox_ref_ECP_pattern_pixels_y.setValue(Ny)


    def _restore_loaded_pattern(self):
        self.measured_ref_ECP = copy.deepcopy(self.measured_ref_ECP_stored)
        Ny, Nx = self.measured_ref_ECP.shape
        self.spinBox_ref_ECP_pattern_pixels_x.setValue(Nx)
        self.spinBox_ref_ECP_pattern_pixels_y.setValue(Ny)
        self.update_display(image=self.measured_ref_ECP, mode='ref_ECCI_measurement')



    # def update_display(self, image, mode='ref_ECCI_measurement'):
    #     image_8bit = image / image.max() * 255
    #     image_to_display = qimage2ndarray.array2qimage(image_8bit.copy())
    #     if mode=='ref_ECCI_measurement':
    #         self.label_image_ref_ECCI_measurement.setPixmap(QtGui.QPixmap(image_to_display))
    #     elif mode=='ref_ECCI_simulation':
    #         self.label_image_ref_ECCI_simulation.setPixmap(QtGui.QPixmap(image_to_display))
    #     elif mode=='Si_ECCI_difference':
    #         self.label_image_Si_ECCI_difference.setPixmap(QtGui.QPixmap(image_to_display))
    #     else:
    #         self.label_messages.setText('No image acquired')







    def _change_angle_step(self):
        number_of_steps = 2 * self.doubleSpinBox_angle_range.value() / float(self.comboBox_angle_step.currentText())
        print(f'number of angle steps = {number_of_steps}')
        self.horizontalScrollBar_tilt_x.setMaximum(int(number_of_steps))
        self.horizontalScrollBar_tilt_y.setMaximum(int(number_of_steps))


    def _set_tilt(self, selector=1, plot=True):
        dict = {1 : {'spinBox' : self.doubleSpinBox_tilt_x,
                     'number_of_steps' : self.horizontalScrollBar_tilt_x,
                     'scrollBar' : self.horizontalScrollBar_tilt_x},
                2: {'spinBox': self.doubleSpinBox_tilt_y,
                    'number_of_steps': self.horizontalScrollBar_tilt_y,
                    'scrollBar': self.horizontalScrollBar_tilt_y}
        }
        if selector==1 or selector==2:
            #if not self.horizontalScrollBar_tilt_x.isSliderDown():
            value = dict[selector]['scrollBar'].value()
            value = utils.map_scrollbar_to_value(x_start=-1*self.doubleSpinBox_angle_range.value(),
                                                 x_end  =+1*self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=dict[selector]['number_of_steps'].maximum(),
                                                 value_to_map=value)
            dict[selector]['spinBox'].setValue(value)

            self.tilt_x = self.doubleSpinBox_tilt_x.value()
            self.tilt_y = self.doubleSpinBox_tilt_y.value()

            if plot: self.calculate_simulated_ECP_pattern()


    def _set_tilt2(self, selector=1, plot=True):
        dict = {1 : {'spinBox' : self.doubleSpinBox_tilt_x,
                     'number_of_steps' : self.horizontalScrollBar_tilt_x,
                     'scrollBar' : self.horizontalScrollBar_tilt_x},
                2: {'spinBox': self.doubleSpinBox_tilt_y,
                    'number_of_steps': self.horizontalScrollBar_tilt_y,
                    'scrollBar': self.horizontalScrollBar_tilt_y}
        }
        if selector==1 or selector==2:
            value = dict[selector]['spinBox'].value()
            value = utils.map_value_to_scrollbar(x_start=-1*self.doubleSpinBox_angle_range.value(),
                                                 x_end  =+1*self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=dict[selector]['scrollBar'].maximum(),
                                                 value_to_map=value)
            dict[selector]['scrollBar'].setValue(int(value))

            self.tilt_x = self.doubleSpinBox_tilt_x.value()
            self.tilt_y = self.doubleSpinBox_tilt_y.value()

            if plot == True:  self.calculate_simulated_ECP_pattern()



    def _set_tilt3(self, selector=1, plot=True):
        dict = {1 : {'spinBox' : self.doubleSpinBox_tilt_x,
                     'number_of_steps' : self.horizontalScrollBar_tilt_x,
                     'scrollBar' : self.horizontalScrollBar_tilt_x},
                2: {'spinBox': self.doubleSpinBox_tilt_y,
                    'number_of_steps': self.horizontalScrollBar_tilt_y,
                    'scrollBar': self.horizontalScrollBar_tilt_y}
        }
        if selector==1 or selector==2:
            value = dict[selector]['scrollBar'].value()
            value = utils.map_scrollbar_to_value(x_start=-1*self.doubleSpinBox_angle_range.value(),
                                                 x_end  =+1*self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=dict[selector]['number_of_steps'].maximum(),
                                                 value_to_map=value)
            dict[selector]['spinBox'].setValue(value)

            self.tilt_x = self.doubleSpinBox_tilt_x.value()
            self.tilt_y = self.doubleSpinBox_tilt_y.value()

            if plot:  self.calculate_simulated_ECP_pattern()







    def _open_tif_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self,
                                                   "QFileDialog.getOpenFileName()",
                                                   "", "TIF files (*.tif);;TIFF files (*.tiff);;All Files (*)",
                                                   options=options)
        if file_name:
            print(file_name)
            if file_name.lower().endswith('.tif') or file_name.lower().endswith('.tiff'):
                data = utils.load_image(file_name)
                success = True
                return (success, file_name, data)
            # other file format, not tiff, for example numpy array data, or txt format
            else:
                try:
                    data = np.loadtxt(file_name)
                    success = True
                    return (success, file_name, data)
                except:
                    return (False, None, None)
        else:
            # no file selected, return nothing
            return (False, None, None)


    def _open_ctf_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self,
                                                   "QFileDialog.getOpenFileName()",
                                                   "", "ctf files (*.ctf);;h5 files (*.h5);;All Files (*)",
                                                   options=options)
        if file_name:
            print(file_name)
            if file_name.lower().endswith('.ctf') or file_name.lower().endswith('.h5'):
                success = True
                return (success, file_name)
        # other file format
        else:
            # no file selected, return nothing
            return (False, None)


    def _open_master_pattern(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self,
                                                   "QFileDialog.getOpenFileName()",
                                                   "", "h5 files (*.h5);;hdf5 files (*.hdf5);;All Files (*)",
                                                   options=options)
        if file_name:
            print(file_name)
            if file_name.lower().endswith('.h5') or file_name.lower().endswith('.hdf5'):
                success = True
                return (success, file_name)
        # other file format, not tiff, for example numpy array data, or txt format
        else:
            # no file selected, return nothing
            return (False, None)


    def _open_ebsp_stack(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self,
                                                   "QFileDialog.getOpenFileName()",
                                                   "", "ebsp files (*.ebsp);;h5 files (*.h5);;All Files (*)",
                                                   options=options)
        if file_name:
            print(file_name)
            if file_name.lower().endswith('.ebsp'):
                success = True
                return (success, file_name)
        # other file format, not tiff, for example numpy array data, or txt format
        else:
            # no file selected, return nothing
            return (False, None)




    def _abort_clicked(self):
        print('------------ abort clicked --------------')
        self._abort_clicked_status = True



def main():
    app = QtWidgets.QApplication([])
    qt_app = GUIMainWindow()
    app.aboutToQuit.connect(qt_app.disconnect)  # cleanup & teardown
    qt_app.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
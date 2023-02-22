import qtdesigner_files.main_gui as gui_main
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtWidgets import QFileDialog
# import qimage2ndarray

# from importlib import reload  # Python 3.4+
#
# from dataclasses import dataclass

import sys, time, os, glob
import numpy as np
import copy

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as _FigureCanvas
from matplotlib.backends.backend_qt5agg import (
    NavigationToolbar2QT as _NavigationToolbar,
)

import h5py
import hyperspy.api as hs
import kikuchipy as kp
import ECP

import utils



class GUIMainWindow(gui_main.Ui_MainWindow, QtWidgets.QMainWindow):
    def __init__(self):
        super(GUIMainWindow, self).__init__()
        self.setupUi(self)

        self.measured_ECP = None
        self.measured_ECP_stored = None
        self.simulated_ECP = None
        self.simulated_ECP_stored = None
        self.difference = None
        self._abort_clicked_status = False
        self.calibrated_tilt_x = 0
        self.calibrated_tilt_y = 0

        self.ecp_reference = ECP.ECP()


        self.tilt_x = 0.0
        self.tilt_y = 0.0

        self.horizontalScrollBar_tilt_y.setValue(50)
        self.horizontalScrollBar_tilt_x.setValue(50)
        self.doubleSpinBox_tilt_y.setValue(0)
        self.doubleSpinBox_tilt_x.setValue(0)
        self.comboBox_angle_step.setCurrentText('0.1')
        self.comboBox_angle_step_automatic.setCurrentText('0.1')
        number_of_steps = 2*self.doubleSpinBox_angle_range.value() / float(self.comboBox_angle_step.currentText())
        self.horizontalScrollBar_tilt_x.setMaximum(int(number_of_steps))
        self.horizontalScrollBar_tilt_y.setMaximum(int(number_of_steps))

        self.setStyleSheet("""QPushButton {
        border: 1px solid lightgray;
        border-radius: 5px;
        background-color: #e3e3e3;
        }""")
        self.setup_connections()
        self.initialise_image_frames()


    def setup_connections(self):
        self.label_messages.setText('Starting up...')
        self.pushButton_abort_automatic_calibration.clicked.connect(lambda: self._abort_clicked())
        #
        self.horizontalScrollBar_tilt_y.valueChanged.connect(lambda: self._set_tilt(angle_slider=2))
        self.horizontalScrollBar_tilt_x.valueChanged.connect(lambda: self._set_tilt(angle_slider=1))
        self.horizontalScrollBar_tilt_y.sliderReleased.connect(lambda: self._set_tilt3(angle_slider=2))
        self.horizontalScrollBar_tilt_x.sliderReleased.connect(lambda: self._set_tilt3(angle_slider=1))
        #
        self.doubleSpinBox_tilt_y.editingFinished.connect(lambda: self._set_tilt2(angle_num=2))
        self.doubleSpinBox_tilt_x.editingFinished.connect(lambda: self._set_tilt2(angle_num=1))
        self.comboBox_angle_step.currentTextChanged.connect(lambda: self._change_angle_step())
        self.doubleSpinBox_angle_range.editingFinished.connect(lambda: self._change_angle_step())
        #
        self.pushButton_open_file_Si_ECCI_measurement.clicked.connect(lambda: self._open_Si_ECCI_measurement_file())
        self.pushButton_crop_meas_ECP_X.clicked.connect(lambda: self._crop_measured_ECP(mode='X'))
        self.pushButton_crop_meas_ECP_Y.clicked.connect(lambda: self._crop_measured_ECP(mode='Y'))
        self.pushButton_crop_meas_ECP_restore.clicked.connect(lambda: self._restore_loaded_pattern())
        self.pushButton_load_ECP_master_pattern.clicked.connect(lambda: self._load_ref_master_pattern())
        self.pushButton_set_ref_ECP_detector.clicked.connect(lambda: self._update_ecp_ref_settings())
        self.pushButton_load_ref_ctf_file.clicked.connect(lambda: self._load_ctf_file())
        self.pushButton_ref_ECP_display.clicked.connect(lambda: self.calculate_simulated_ECP_pattern())
        self.pushButton_run_automatic_calibration.clicked.connect(lambda: self.run_automatic_calibration())


    def initialise_image_frames(self):
        self.figure_ECCI_exp = plt.figure(10)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_ECCI_exp = _FigureCanvas(self.figure_ECCI_exp)
        self.toolbar_ECCI_exp = _NavigationToolbar(self.canvas_ECCI_exp, self)
        #
        self.label_image_Si_ECCI_measurement.setLayout(QtWidgets.QVBoxLayout())
        self.label_image_Si_ECCI_measurement.layout().addWidget(self.toolbar_ECCI_exp)
        self.label_image_Si_ECCI_measurement.layout().addWidget(self.canvas_ECCI_exp)

        self.figures = {1 : {'fig' : self.figure_ECCI_exp, 'canvas': self.canvas_ECCI_exp, 'toolbar': self.toolbar_ECCI_exp}  }
        ################################################################################################
        self.figure_ECCI_sim = plt.figure(11)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_ECCI_sim = _FigureCanvas(self.figure_ECCI_sim)
        self.toolbar_ECCI_sim = _NavigationToolbar(self.canvas_ECCI_sim, self)
        #
        self.label_image_Si_ECCI_simulation.setLayout(QtWidgets.QVBoxLayout())
        self.label_image_Si_ECCI_simulation.layout().addWidget(self.toolbar_ECCI_sim)
        self.label_image_Si_ECCI_simulation.layout().addWidget(self.canvas_ECCI_sim)

        self.figures[2] = {'fig' : self.figure_ECCI_sim, 'canvas': self.canvas_ECCI_sim, 'toolbar': self.toolbar_ECCI_sim}

        ################################################################################################
        self.figure_ECCI_diff = plt.figure(12)
        plt.axis("off")
        plt.tight_layout()
        plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.01)
        self.canvas_ECCI_diff = _FigureCanvas(self.figure_ECCI_diff)
        self.toolbar_ECCI_diff = _NavigationToolbar(self.canvas_ECCI_diff, self)
        #
        self.label_image_Si_ECCI_difference.setLayout(QtWidgets.QVBoxLayout())
        self.label_image_Si_ECCI_difference.layout().addWidget(self.toolbar_ECCI_diff)
        self.label_image_Si_ECCI_difference.layout().addWidget(self.canvas_ECCI_diff)

        self.figures[3] = {'fig' : self.figure_ECCI_diff, 'canvas': self.canvas_ECCI_diff, 'toolbar': self.toolbar_ECCI_diff}
        ################################################################################################


    def _update_ecp_ref_settings(self):
        if self.measured_ECP is not None:
            energy = self.doubleSpinBox_ref_ECP_energy.value()
            projection = self.comboBox_ref_ECP_projection.currentText()
            hemispheres = self.comboBox_ref_ECP_hemispheres.currentText()
            shape = self.measured_ECP.shape
            pc_x = self.doubleSpinBox_ref_ECP_pc_x.value()
            pc_y = self.doubleSpinBox_ref_ECP_pc_y.value()
            pc_z = self.doubleSpinBox_ref_ECP_pc_z.value()
            pixel_size = self.doubleSpinBox_ref_ECP_pixel_size.value()
            tilt = self.doubleSpinBox_ref_ECP_detector_tilt.value()
            sample_tilt = self.doubleSpinBox_ref_ECP_sample_tilt.value()
            binning = self.spinBox_ref_ECP_binning.value()
            self.ecp_reference.update_settings(self,
                                               energy=energy,
                                               pc_x = pc_x,
                                               pc_y = pc_y,
                                               pc_z = pc_z,
                                               pixel_size = pixel_size,
                                               binning = binning,
                                               tilt = tilt,
                                               sample_tilt = sample_tilt,
                                               projection = projection,
                                               hemispheres = hemispheres,
                                               detector_shape = shape)
            self.label_messages.setText('ref ECP detector settings uploaded')
            print(energy, pc_x, pc_y, pc_z, pixel_size, binning, tilt, sample_tilt, projection, hemispheres, shape)
        else:
            self.label_messages.setText('First, load a measured ECP/ECCI pattern to define the number of pixels in the pattern')



    def _open_Si_ECCI_measurement_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self,
                                                   "QFileDialog.getOpenFileName()",
                                                   "", "TIF files (*.tif);;TIFF files (*.tiff);;All Files (*)",
                                                   options=options)
        if file_name:
            print(file_name)
            self.label_messages.setText('ref ECP measurement filename: ' + file_name)
            if file_name.lower().endswith('.tif') or file_name.lower().endswith('.tiff'):
                self.measured_ECP = utils.load_image(file_name)
                Ny, Nx = self.measured_ECP.shape
                self.spinBox_ref_ECP_pattern_pixels_x.setValue(Nx)
                self.spinBox_ref_ECP_pattern_pixels_y.setValue(Ny)
                self.update_display(image=self.measured_ECP, mode='Si_ECCI_measurement')

            # other file format, not tiff, for example numpy array data, or txt format
            else:
                try:
                    self.measured_ECP = np.loadtxt(file_name)
                    self.update_display(image=self.measured_ECP, mode='Si_ECCI_measurement')
                except:
                    self.label_messages.setText('File or mode not supported')


    def _load_ref_master_pattern(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self,
                                                   "QFileDialog.getOpenFileName()",
                                                   "", "h5 files (*.h5);;hdf5 files (*.hdf5);;All Files (*)",
                                                   options=options)
        if file_name:
            print(file_name)
            self.label_messages.setText('ref ECP measurement filename: ' + file_name)
            self.label_ecp_master_pattern_path.setText(file_name)
            self._update_ecp_ref_settings()
            if file_name.lower().endswith('.h5') or file_name.lower().endswith('.hdf5'):
                output = \
                    self.ecp_reference.load_master_pattern(path_to_master_pattern=file_name)
                self.label_messages.setText(output)


    def _load_ctf_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self,
                                                   "QFileDialog.getOpenFileName()",
                                                   "", "ctf files (*.ctf);;h5 files (*.h5);;All Files (*)",
                                                   options=options)
        if file_name:
            print(file_name)
            self.label_messages.setText('ctf filename: ' + file_name)
            self.label_ref_ctf_file.setText(file_name)
            self._update_ecp_ref_settings()
            if file_name.lower().endswith('.ctf'):
                output = \
                    self.ecp_reference.load_xmap(file_name=file_name)
                self.label_messages.setText(str(output))


    def calculate_simulated_ECP_pattern(self, plot=True):
        try:
            self.simulated_ECP = \
                self.ecp_reference.calculate_ecp_pattern(tilt_x=self.tilt_x,
                                                         tilt_y=self.tilt_y)
            self.update_display(image=self.simulated_ECP,
                                mode='Si_ECCI_simulation')

            # if both experimental pattern and similated pattern exist,
            # calculate their difference and plot if True
            if (self.simulated_ECP is not None) and (self.measured_ECP is not None):
                self.difference = utils.calculate_difference(image1=self.measured_ECP,
                                                             image2=self.simulated_ECP)
                #calculate the similarity of two images, display MAC
                self.mac = utils.modal_assurance_criterion(image1=self.measured_ECP,
                                                           image2=self.simulated_ECP)
                self.doubleSpinBox_mac.setValue(self.mac)
                if plot==True:
                    self.update_display(image=self.difference,
                                        mode='difference')
        except Exception as e:
            print('could not simulate the pattern', e)


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
                self._set_tilt2(angle_num=1, plot=False) # update the slider
                self._set_tilt2(angle_num=2, plot=False) # update the slider
                _plot_ = self.checkBox_plot_while_running.isChecked()
                self.calculate_simulated_ECP_pattern(plot=_plot_)
                if self.mac > mac_max:
                    mac_max = self.mac
                    tilt_x_max = tilt_x
                    tilt_y_max = tilt_y_max
                    self.doubleSpinBox_tilt_x_calibrated.setValue(tilt_x_max)
                    self.doubleSpinBox_tilt_y_calibrated.setValue(tilt_y_max)
                self.repaint()
                QtWidgets.QApplication.processEvents()

        # update the calibrated tilt values at the end of the run
        self.calibrated_tilt_x = tilt_x_max
        self.calibrated_tilt_y = tilt_y_max
        self.label_messages.setText('Calibration completed')





    def _crop_measured_ECP(self, mode='X'):
        if self.measured_ECP is not None:
            start = self.spinBox_ref_meas_ECP_crop_start.value()
            end = self.spinBox_ref_meas_ECP_crop_end.value()
            self.measured_ECP_stored = copy.deepcopy(self.measured_ECP)
            if mode=="X":
                self.measured_ECP = self.measured_ECP[:, start : end]
            else:
                self.measured_ECP = self.measured_ECP[start : end, :]
            self.update_display(image=self.measured_ECP, mode='Si_ECCI_measurement')
            Ny, Nx = self.measured_ECP.shape
            self.spinBox_ref_ECP_pattern_pixels_x.setValue(Nx)
            self.spinBox_ref_ECP_pattern_pixels_y.setValue(Ny)


    def _restore_loaded_pattern(self):
        self.measured_ECP = copy.deepcopy(self.measured_ECP_stored)
        Ny, Nx = self.measured_ECP.shape
        self.spinBox_ref_ECP_pattern_pixels_x.setValue(Nx)
        self.spinBox_ref_ECP_pattern_pixels_y.setValue(Ny)
        self.update_display(image=self.measured_ECP, mode='Si_ECCI_measurement')



    # def update_display(self, image, mode='Si_ECCI_measurement'):
    #     image_8bit = image / image.max() * 255
    #     image_to_display = qimage2ndarray.array2qimage(image_8bit.copy())
    #     if mode=='Si_ECCI_measurement':
    #         self.label_image_Si_ECCI_measurement.setPixmap(QtGui.QPixmap(image_to_display))
    #     elif mode=='Si_ECCI_simulation':
    #         self.label_image_Si_ECCI_simulation.setPixmap(QtGui.QPixmap(image_to_display))
    #     elif mode=='Si_ECCI_difference':
    #         self.label_image_Si_ECCI_difference.setPixmap(QtGui.QPixmap(image_to_display))
    #     else:
    #         self.label_messages.setText('No image acquired')


    def update_display(self, image, mode='Si_ECCI_measurement'):
        cmap = 'gray'
        if mode=='Si_ECCI_measurement':
            key = 1
        if mode=='Si_ECCI_simulation':
            key = 2
        if mode=='difference':
            key = 3
            cmap = 'bwr'

        self.figures[key]['fig'].clear()
        self.figures[key]['fig'].patch.set_facecolor(
            (240 / 255, 240 / 255, 240 / 255))
        self.ax = self.figures[key]['fig'].add_subplot(111)
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)
        self.ax.imshow(image, cmap=cmap)
        self.figures[key]['canvas'].draw()





    def _change_angle_step(self):
        number_of_steps = 2 * self.doubleSpinBox_angle_range.value() / float(self.comboBox_angle_step.currentText())
        print(f'number of angle steps = {number_of_steps}')
        self.horizontalScrollBar_tilt_x.setMaximum(int(number_of_steps))
        self.horizontalScrollBar_tilt_y.setMaximum(int(number_of_steps))


    def _set_tilt(self, angle_slider=1):
        if angle_slider==1:
            if not self.horizontalScrollBar_tilt_x.isSliderDown():
                pass #print("updating")
            value = self.horizontalScrollBar_tilt_x.value()
            value = utils.map_scrollbar_to_value(x_start=-1*self.doubleSpinBox_angle_range.value(), x_end=+1*self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=self.horizontalScrollBar_tilt_x.maximum(),
                                                 value_to_map=value)
            self.doubleSpinBox_tilt_x.setValue(value)
            self.tilt_x = value
            #self.calculate_simulated_ECP_pattern()
        elif angle_slider==2:
            if not self.horizontalScrollBar_tilt_y.isSliderDown():
                pass #print("i am updating")
            value = self.horizontalScrollBar_tilt_y.value()
            value = utils.map_scrollbar_to_value(x_start=-1 * self.doubleSpinBox_angle_range.value(),
                                                 x_end=+1 * self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=self.horizontalScrollBar_tilt_y.maximum(),
                                                 value_to_map=value)
            self.doubleSpinBox_tilt_y.setValue(value)
            self.tilt_y = value
            #self.calculate_simulated_ECP_pattern()


    def _set_tilt2(self, angle_num=1, plot=True):
        if angle_num==1:
            value = self.doubleSpinBox_tilt_x.value()
            self.tilt_x = value
            value = utils.map_value_to_scrollbar(x_start=-1*self.doubleSpinBox_angle_range.value(),
                                                 x_end=+1*self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=self.horizontalScrollBar_tilt_x.maximum(),
                                                 value_to_map=value)
            self.horizontalScrollBar_tilt_x.setValue(int(value))
            self.calculate_simulated_ECP_pattern()
        elif angle_num==2:
            value = self.doubleSpinBox_tilt_y.value()
            self.tilt_y = value
            value = utils.map_value_to_scrollbar(x_start=-1*self.doubleSpinBox_angle_range.value(),
                                                 x_end=+1*self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=self.horizontalScrollBar_tilt_y.maximum(),
                                                 value_to_map=value)
            self.horizontalScrollBar_tilt_y.setValue(int(value))
        if plot==True:
            self.calculate_simulated_ECP_pattern()


    def _set_tilt3(self, angle_slider=1, plot=True):
        if angle_slider==1:
            value = self.horizontalScrollBar_tilt_x.value()
            value = utils.map_scrollbar_to_value(x_start=-1*self.doubleSpinBox_angle_range.value(), x_end=+1*self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=self.horizontalScrollBar_tilt_x.maximum(),
                                                 value_to_map=value)
            self.doubleSpinBox_tilt_x.setValue(value)
            self.tilt_x = value
        elif angle_slider==2:
            value = self.horizontalScrollBar_tilt_y.value()
            value = utils.map_scrollbar_to_value(x_start=-1 * self.doubleSpinBox_angle_range.value(),
                                                 x_end=+1 * self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=self.horizontalScrollBar_tilt_y.maximum(),
                                                 value_to_map=value)
            self.doubleSpinBox_tilt_y.setValue(value)
            self.tilt_y = value
        if plot:
            self.calculate_simulated_ECP_pattern()


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
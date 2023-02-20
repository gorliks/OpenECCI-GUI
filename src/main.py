import qtdesigner_files.main_gui as gui_main
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtWidgets import QFileDialog
# import qimage2ndarray

# from importlib import reload  # Python 3.4+
#
# from dataclasses import dataclass

import sys, time, os, glob
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as _FigureCanvas
from matplotlib.backends.backend_qt5agg import (
    NavigationToolbar2QT as _NavigationToolbar,
)

import h5py
import hyperspy.api as hs
# import kikuchipy as kp

import utils



class GUIMainWindow(gui_main.Ui_MainWindow, QtWidgets.QMainWindow):
    def __init__(self):
        super(GUIMainWindow, self).__init__()
        self.setupUi(self)

        self.calibration_tilt_x = 0.0
        self.calibration_tilt_y = 0.0

        self.horizontalScrollBar_tilt_y.setValue(50)
        self.horizontalScrollBar_tilt_x.setValue(50)
        self.doubleSpinBox_tilt_y.setValue(0)
        self.doubleSpinBox_tilt_x.setValue(0)
        self.comboBox_angle_step.setCurrentText('0.1')
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



    def _open_Si_ECCI_measurement_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self,
                                                   "QFileDialog.getOpenFileName()",
                                                   "", "TIF files (*.tif);;TIFF files (*.tiff);;All Files (*)",
                                                   options=options)
        if file_name:
            print(file_name)
            if file_name.lower().endswith('.tif') or file_name.lower().endswith('.tiff'):
                self.image = utils.load_image(file_name)
                print('self.image: ', self.image.shape)
                self.update_display(image=self.image, mode='Si_ECCI_measurement')

            # other file format, not tiff, for example numpy array data, or txt format
            else:
                try:
                    self.image = np.loadtxt(file_name)
                    self.update_image(image=self.image)
                except:
                    self.label_messages.setText('File or mode not supported')



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
        if mode=='Si_ECCI_measurement':
            key = 1

        self.figures[key]['fig'].clear()
        self.figures[key]['fig'].patch.set_facecolor(
            (240 / 255, 240 / 255, 240 / 255))
        self.ax = self.figures[key]['fig'].add_subplot(111)
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)
        self.ax.imshow(image, cmap='gray')
        self.figures[key]['canvas'].draw()





    def _change_angle_step(self):
        number_of_steps = 2 * self.doubleSpinBox_angle_range.value() / float(self.comboBox_angle_step.currentText())
        print(f'number of angle steps = {number_of_steps}')
        self.horizontalScrollBar_tilt_x.setMaximum(int(number_of_steps))
        self.horizontalScrollBar_tilt_y.setMaximum(int(number_of_steps))

    def _set_tilt(self, angle_slider=1):
        if angle_slider==1:
            if self.horizontalScrollBar_tilt_x.isSliderDown():
                pass
            else:
                print("updating")
            value = self.horizontalScrollBar_tilt_x.value()
            value = utils.map_scrollbar_to_value(x_start=-1*self.doubleSpinBox_angle_range.value(), x_end=+1*self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=self.horizontalScrollBar_tilt_x.maximum(),
                                                 value_to_map=value)
            self.doubleSpinBox_tilt_x.setValue(value)
        elif angle_slider==2:
            if self.horizontalScrollBar_tilt_y.isSliderDown():
                pass
            else:
                print("i am updating")
            value = self.horizontalScrollBar_tilt_y.value()
            value = utils.map_scrollbar_to_value(x_start=-1 * self.doubleSpinBox_angle_range.value(),
                                                 x_end=+1 * self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=self.horizontalScrollBar_tilt_y.maximum(),
                                                 value_to_map=value)
            self.doubleSpinBox_tilt_y.setValue(value)

    def _set_tilt2(self, angle_num=1):
        if angle_num==1:
            value = self.doubleSpinBox_tilt_x.value()
            value = utils.map_value_to_scrollbar(x_start=-1*self.doubleSpinBox_angle_range.value(),
                                                 x_end=+1*self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=self.horizontalScrollBar_tilt_x.maximum(),
                                                 value_to_map=value)
            self.horizontalScrollBar_tilt_x.setValue(int(value))
        elif angle_num==2:
            value = self.doubleSpinBox_tilt_y.value()
            value = utils.map_value_to_scrollbar(x_start=-1*self.doubleSpinBox_angle_range.value(),
                                                 x_end=+1*self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=self.horizontalScrollBar_tilt_y.maximum(),
                                                 value_to_map=value)
            self.horizontalScrollBar_tilt_y.setValue(int(value))

    def _set_tilt3(self, angle_slider=1):
        if angle_slider==1:
            value = self.horizontalScrollBar_tilt_x.value()
            value = utils.map_scrollbar_to_value(x_start=-1*self.doubleSpinBox_angle_range.value(), x_end=+1*self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=self.horizontalScrollBar_tilt_x.maximum(),
                                                 value_to_map=value)
        elif angle_slider==2:
            value = self.horizontalScrollBar_tilt_y.value()
            value = utils.map_scrollbar_to_value(x_start=-1 * self.doubleSpinBox_angle_range.value(),
                                                 x_end=+1 * self.doubleSpinBox_angle_range.value(),
                                                 number_of_steps=self.horizontalScrollBar_tilt_y.maximum(),
                                                 value_to_map=value)






def main():
    app = QtWidgets.QApplication([])
    qt_app = GUIMainWindow()
    app.aboutToQuit.connect(qt_app.disconnect)  # cleanup & teardown
    qt_app.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

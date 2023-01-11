import serial
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore
from interpolation import inverse

g = 6
ser = serial.Serial("/dev/cu.usbserial-AR0KL1IW",2000000)
app = pg.mkQApp("Vision thermique")
window = pg.GraphicsLayoutWidget()
view = window.addViewBox()
img = pg.ImageItem()
timer = QtCore.QTimer()

window.show()
window.setWindowTitle("Vision thermique interpolée par inverse des distances")
view.setAspectLocked(True)
view.addItem(img)
view.setRange(QtCore.QRectF(0, 0, 8*g, 8*g))
timer.setSingleShot(True)

def mat():
    line = ser.readline().decode("utf-8")
    matrice = np.fromstring(line,dtype=float,sep=" ").reshape(8,8)
    interpo = np.array(inverse(matrice,g))
    return interpo

def update():
    img.setImage(mat())
    timer.start(1)

timer.timeout.connect(update)
update()

if __name__ == "__main__":
    pg.exec()

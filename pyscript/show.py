from matplotlib import pyplot as plt
import struct
import numpy as np
def r1():
    fURL = "qumian.txt"
    m = np.loadtxt(fURL)
    plt.imshow(m)
    plt.colorbar()
    plt.show()
if __name__ == "__main__":
    src_URL = "..//models//RT_SynMod2D.src"
    rec_FW_URL = "..//models//RT_SynMod2D_FW.rec"
    rec_INV_URL = "..//models//RT_SynMod2D_INV.rec"
    ref_URL = "..//models//RT_SynMod2D.ref"
    vel_FW_URL = "..//models//RT_SynMod2D_FW.vp"

    weight_INV_URL = "..//models//RT_SynMod2D_FW.vp"

    src_test = "..//RT_SynMod2D.src"
    vel_FW_URL = "..//mod//RT_SynMod2D_FW.vp"
    showURL = "..//mod//RT_SynMod2D_FW_Src9_TTd.2d"
    
    
    fURL = "..//mod//RT_SynMod2D.finalvp"
    vel_INV_URL = "..//mod//RT_SynMod2D_INV.vp"

    with open(fURL, "rb") as fstream:
        k = fstream.read( )
    print(len(k))

    dpl = 8 # double precision length
    ylen = 500 # 160
    xlen = 500 # 600
    resMap = []
    for i in range(int(len(k)/dpl)):
        start = i * 8
        end = (i+1) * 8
        val = (struct.unpack('d', k[start:end])) # little-end
        resMap.append(val)
    showMap = np.array(resMap).reshape([ylen,xlen])
    plt.imshow(showMap)
    plt.colorbar()
    plt.show()
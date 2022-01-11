import struct
import numpy
from matplotlib import pyplot as plt

if __name__ == "__main__":
    # Forward model 
    vel_FW_URL = "..//models//RT_SynMod2D_FW.vp"
    src_URL = "..//models//RT_SynMod2D.src"
    rec_FW_URL = "..//models//RT_SynMod2D_FW.rec"
    ref_URL = "..//models//RT_SynMod2D.ref"
    #wr_
    # Inversion model
    rec_INV_URL = "..//models//RT_SynMod2D_Inv.wr"
    
    # Test Model
    showURL2 = "..//mod/RT_SynMod2D_FW_Tr2D.bin"
    src_test = "..//RT_SynMod2D.src"
    fURL = rec_INV_URL

    with open(fURL, "rb") as fstream:
        k = fstream.read( )
    print(len(k)/8)

    for i in range(8*1):
        start = i * 8
        end = (i+1) * 8
        #res = struct.unpack('d', k[start:end]
        #print(struct.unpack('d', k[start:end]))
    for i in range(int(len(k)/8)):
        start = i * 8
        end = (i+1) * 8
        res = struct.unpack('d', k[start:end])
        if not res[0] == 0.5:
            print("x")
    #print(struct.unpack('d', k[start:end]))
    # print(k)
    # 
    print("A.O.H.")

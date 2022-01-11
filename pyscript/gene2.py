import numpy as np
import struct


if __name__ == "__main__":
    rec_FW_URL = "..//mod//RT_SynMod2D_FW.rec"
    rec_TIME_URL = "..//mod//RT_SynMod2D_FW_Tr2D.bin"

    rec_INV_URL = "..//mod//RT_SynMod2D_INV.rec"
    vel_INV_URL = "..//mod//RT_SynMod2D_INV.vp"
    weight_INV_URL = "..//mod//RT_SynMod2D_Inv.wr"

    with open(rec_FW_URL, "rb") as fstream:
        k = fstream.read( )
    fstream.close()
    print(len(k)/8)
    data = []
    for i in range(int(len(k)/8)): 
        start = i * 8
        end = (i+1) * 8
        data.append(struct.unpack('d', k[start:end]))
    totalTrace = int(len(k)/8/4)
    dataList = np.array(data).reshape([totalTrace,4]) # 4xn 行矩阵
    print(dataList[0])
    
    # ----------------- open time file
    with open(rec_TIME_URL, "rb") as fstream:
        k = fstream.read( )
    fstream.close()
    data = []
    for i in range(int(len(k)/8)): 
        start = i * 8
        end = (i+1) * 8
        data.append(struct.unpack('d', k[start:end]))
    print(len(data[0]),len(dataList[0]))
    res = []
    for i in range(len(dataList)):
        res.append([dataList[i][0], dataList[i][1],dataList[i][2],dataList[i][3], data[i*2][0], data[i*2+1][0]])
    
    # ------------------ WRITE time into record
    fstream = open(rec_INV_URL, "wb")
    for i in range(len(res)):
        for j in range(len(res[0])):
            temp = res[i][j]
            fstream.write(struct.pack('d', temp))
    fstream.close()
    print(res[0])
    
    # ------------------ weight value
    fstream = open(weight_INV_URL, "wb")
    for i in range(500):
        for j in range(500):
            temp = 0.5 # if c1 choose 1
            fstream.write(struct.pack('d', temp))
    fstream.close()

    # ------------------ velocity ini model
    fstream = open(vel_INV_URL, "wb")
    for i in range(500):
        for j in range(500):
            temp = 1500 + i * 5 * 1
            fstream.write(struct.pack('d', temp))
    fstream.close()
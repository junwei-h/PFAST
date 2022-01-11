import numpy as np
import struct

def writeFile2D(sm, fURL):
    fstream = open(fURL, 'wb')
    for i in range(len(m)):
        for j in range(len(m[0])):
            res = m[i][j]
            fstream.write(struct.pack('d', res))
    fstream.close()
    return 0

if __name__ == "__main__":

    # -------------- vp file
    fURL = "qumian.txt"
    m = np.loadtxt(fURL)    
    #struct.unpack('d', k[start:end])) # little-end
    vpURL = "./mod/RT_SynMod2D_FW.vp"
    print(len(m), len(m[0]))
    fstream = open(vpURL, 'wb')
    for i in range(len(m)):
        for j in range(len(m[0])):
            res = m[i][j]
            fstream.write(struct.pack('d', res))
    fstream.close()

    fstream = open(vpURL, 'wb')
    for i in range(500):
        for j in range(500):
            #if i * 5 < 1000:
            #    res = 1540 - i * 5 * 0.04 
            #else:
            #    if i * 5 < 1100:
            #        res = 1500
            #    else:
            #        res = 1500 + i * 5 * 0.04
            res = 1500 + i * 5 * 1
            if 200< i <300:
                if 200 < j < 300:
                    res = 2000

            fstream.write(struct.pack('d', res))
    fstream.close()

    # -------------- source file
    m = []
    src_x = 40 # srcNum
    src_x_interval = 50 # add grid nums
    traceNum = 400
    y = 5 # do not set to 0 due to the height

    for i in range(src_x):
        x = i * src_x_interval + 400 # from 0 to 2450
        m.append([y,x,traceNum])

    srcURL = "./mod/RT_SynMod2D.src"
    fstream = open(srcURL, 'wb')
    for i in range(len(m)):
        for j in range(len(m[0])):
            res = m[i][j]
            fstream.write(struct.pack('d', res))
    fstream.close()

    # -------------- recorder file
    recURL = "./mod/RT_SynMod2D_FW.rec"
    binURL = "./mod/RT_SynMod2D_FW.rec"
    m = []

    rec_x_interval = 5
    start_rec = 200
    normal_y = -1
    normal_x = 0
    y = 5 #
    for j in range(src_x): # shot iteration
        for i in range(traceNum):
            x = i*rec_x_interval + start_rec  # from 0 to 490
            m.append([y,x,normal_y, normal_x])

    fstream = open(recURL, 'wb')
    for i in range(len(m)):
        for j in range(len(m[0])):
            res = m[i][j]
            fstream.write(struct.pack('d', res))
    fstream.close()

    #-------------- reflector file
    refURL = "./mod/RT_SynMod2D.ref"
    m = []
    ref_x = 500
    for i in range(ref_x): # r1
        m.append([2500,-1,0])
    #for i in range(ref_x): # r1
    #    m.append([2400,-1,0])
    #for i in range(ref_x/2): # r2
        #m.append([500,-1,0])

    fstream = open(refURL, 'wb')
    for i in range(len(m)):
        for j in range(len(m[0])):
            res = m[i][j]
            fstream.write(struct.pack('d', res))
    fstream.close()   
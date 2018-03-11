#!/bin/bash

mpiCC -D_FILE_OFFSET_BITS=64 -DNoInverse RT_PFAST2D.cpp fftn.cpp -o ../bin/RT_PFAST2D_FW
mpiCC -D_FILE_OFFSET_BITS=64 -DNoInverse -DNoShotGather RT_PFAST2D.cpp fftn.cpp -o ../bin/RT_PFAST2D_FW_NoS

mpiCC -D_FILE_OFFSET_BITS=64 RT_PFAST2D.cpp fftn.cpp -o ../bin/RT_PFAST2D
mpiCC -D_FILE_OFFSET_BITS=64 RT_PFAST2D_BENCH.cpp fftn.cpp -o ../bin/RT_PFAST2D_BENCH

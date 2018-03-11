#!/bin/bash

mpiCC -D_FILE_OFFSET_BITS=64 -DNoInverse T_PFAST3D.cpp fftn.cpp -o ../bin/T_PFAST3D_FW

mpiCC -D_FILE_OFFSET_BITS=64 T_PFAST3D.cpp fftn.cpp -o ../bin/T_PFAST3D
mpiCC -D_FILE_OFFSET_BITS=64 T_PFAST3D_BENCH.cpp fftn.cpp -o ../bin/T_PFAST3D_BENCH


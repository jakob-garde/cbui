#!/bin/sh
cd "$(dirname "$0")"
cd lib
ld -r -b binary -o all_res.o all.res
cd ..
g++ main.cpp -o cbui -lGL -lGLEW -lglfw lib/all_res.o
g++ -g main.cpp -o cbui_dbg -lGL -lGLEW -lglfw lib/all_res.o

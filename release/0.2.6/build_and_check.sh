#!/bin/sh
ld -r -b binary -o all_res.o all.res
g++ main.cpp -o cbui_version_check -lGL -lGLEW -lglfw all_res.o
rm all_res.o
./cbui_version_check

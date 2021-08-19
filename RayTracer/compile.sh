#!/bin/sh
g++ -O4 -g svdDynamic.c RayTracer.c -fopenmp utils.c -lm -o RayTracer
if %errorlevel% != 0 echo "Did not Compiled"
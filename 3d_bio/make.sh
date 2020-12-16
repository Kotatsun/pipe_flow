#!/bin/bash
ifort -qopenmp -mcmodel=large -shared-intel $1.f90 -o $1.out
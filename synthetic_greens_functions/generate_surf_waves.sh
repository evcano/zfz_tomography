#!/bin/bash

# source depth
HS=1.2
# receiver depth
HR=1.2
# number of modes
NMOD=4

# clean files
rm *.lov *.ray *.PLT *.egn *.out *.dat
rm ./all_modes/*
rm ./mode_0/*
rm ./mode_1/*
rm ./mode_2/*

# generate synthetics
sprep96 -M MODEL.TXT -d DIST.TXT -HS $HS -HR $HR -NMOD $NMOD -L -R

sdisp96 > sdisp96.out

slegn96 -NOQ > slegn96.out
sregn96 -NOQ > sregn96.out


spulse96 -d DIST.TXT -D -i > all_modes.out
spulse96 -d DIST.TXT -D -i -FUND > mode_0.out
spulse96 -d DIST.TXT -D -i -M 1 > mode_1.out
spulse96 -d DIST.TXT -D -i -M 2 > mode_2.out

# generate SAC files
cd ./all_modes
f96tosac -FMT 1 < ../all_modes.out
cd ..

cd ./mode_0
f96tosac -FMT 1 < ../mode_0.out
cd ..

cd ./mode_1
f96tosac -FMT 1 < ../mode_1.out
cd ..

cd ./mode_2
f96tosac -FMT 1 < ../mode_2.out
cd ..

# plot phase vel dispersion curves
sdpsrf96 -L -ASC
sdpsrf96 -R -ASC

# plot group vel dispersion curves
sdpegn96 -L -U -ASC
sdpegn96 -R -U -ASC

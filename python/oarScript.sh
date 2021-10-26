#!/bin/bash

#--- set parameters
lmpDataFile=./dataPeriodicATOM.txt
lx=80.0
ly=80
hards=-1.9
fyild=-1.0e-04
K=2.0
Gbar=1.0
pcrit=1.0e+06
rsidp=0.0
redup=0.0
FRICT=65.0
reduf=0.0
rsidf=1.0
mnunx=0.1
reduc=0.0
resid=0.1
pinit=8.0
UNFRM=1
ZSHER=0
crit='von-mises'
shear='simple'
rcut=4.0

python init.py $lmpDataFile $lx $ly $hards $fyild $K $Gbar $pcrit $rsidp $redup $FRICT $reduf $rsidf $mnunx $reduc $resid $pinit $UNFRM $ZSHER $crit $shear $rcut 


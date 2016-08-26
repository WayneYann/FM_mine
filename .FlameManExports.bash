### change the location of the FlameMaster package here
export FlameManPath=~/FlameMaster

### do not change the following couple of lines
export FlameManSource=$FlameManPath/FlameManSource
export FlameManLib=$FlameManPath/FlameManLibs
export myData=$FlameManPath/FlameManData
alias FlameMaster=$FlameManLib/FlameMaster
alias LT='$FlameManLib/ListTool -M'
export PATH=$PATH:$FlameManLib

### You have to choose between two possible options for FLEXPATH depending
### on weather flex is installed on your machine and in the current path or not.
### You can find out by typing 'which flex' in the command line. If you get a valid
### path as a result, this means flex is installed. Then choose the first option.
### If flex cannot be found, choose the second one. Then, flex will be installed
### during the FlameMaster installation. You will do this following the procedure
### outlined in the Readme file.
###
### Option 1: flex already exists on your machine (don't insert the path in the FLEXPATH definition)
export FLEXPATH=""
### Option 2: flex will be installed to $FlameManLib
#export FLEXPATH=$FlameManLib/bin/


### choose one of the following sets

export FlameManMach=LINUXGXX
export FlameManFF=ifort
export FlameManF90=ifort
export FlameManCC=icc
export FlameManCCOpts="-O3 -xAVX -ip"
export FlameManCXX=icpc
export FlameManCXXOpts="-O3 -xAVX -ip"
export FlameManF77Opts="-O3 -xAVX -ip"
export FlameManLIBF77=""
export FlameManLDFortOpts="-lifcore -limf -lsvml -lm -lipgo -lirc -lirc_s -ldl"
export FlameManYACC="yacc"
export FlameManLDCCOpts=""
export FlameManRANLIB='ranlib $(LIBRARY)'
export FlameManCVODEDIR="/opt/sundials/include"
export FlameManCVODEVERSION="SUNDIALS24"

# following for Mac OSX
#export FlameManMach=DIGITAL
#export FlameManFF=fc
#export FlameManCC=cc
#export FlameManCCOpts="-O3"
#export FlameManCXX=CC
#export FlameManCXXOpts="-O3"
#export FlameManF77Opts="-O3"
#export FlameManLIBF77="/sw/lib/libg2c.a"
#export FlameManLDFortOpts=""
#export FlameManLDCCOpts="-lm -lstdc++"
#export FlameManRANLIB='ranlib $(LIBRARY)'

#export FlameManMach=DIGITAL
#export FlameManFF=f77
#export FlameManCC=cc
#export FlameManCCOpts=""
#export FlameManCXX=CC
#export FlameManCXXOpts=""
#export FlameManF77Opts=""
#export FlameManLIBF77=""
#export FlameManLDFortOpts=""
#export FlameManLDCCOpts="-lm -lUfor -lfor -lFutil -lots"
#export FlameManRANLIB=""

#export FlameManMach=SUNGXX
#export FlameManFF=f77
#export FlameManCC=gcc
#export FlameManCCOpts="-ansi"
#export FlameManCXX=g++
#export FlameManCXXOpts="-x c++"
#export FlameManF77Opts=""
#export FlameManLIBF77="/optx/SUNWspro/SC4.0/lib/libsunmath.a /optx/SUNWspro/SC4.0/lib/libF77.a /optx/SUNWspro/SC4.0/lib/libm.a"
#export FlameManLDFortOpts=""
#export FlameManLDCCOpts=""
#export FlameManRANLIB='ranlib $(LIBRARY)'

#export FlameManMach=SUN
#export FlameManFF=f77
#export FlameManCC=cc
#export FlameManCCOpts="-Xc"
#export FlameManCXX=CC
#export FlameManCXXOpts=""
#export FlameManF77Opts=""
#export FlameManLIBF77="/usr/lang/SC1.0/libF77.a /usr/lang/SC1.0/libm.a"
#export FlameManLDFortOpts=""
#export FlameManLDCCOpts=""
#export FlameManRANLIB=""

# choose following for SGI
#export FlameManMach=SUN
#export FlameManFF=f77
#export FlameManCC=cc
#export FlameManCCOpts="-n32"
#export FlameManCXX=CC
#export FlameManCXXOpts="-n32"
#export FlameManF77Opts="-n32"
#export FlameManLIBF77=""
#export FlameManLDFortOpts="-lftn"
#export FlameManLDCCOpts="-n32"
#export FlameManRANLIB=""

#export FlameManMach=HP
#export FlameManFF=f77
#export FlameManCC=cc
#export FlameManCCOpts="-Aa -c -DHP -z"
#export FlameManCXX=CC
#export FlameManCXXOpts="-Aa -DHP -z +a1"
#export FlameManF77Opts="-U"
#export FlameManLIBF77=""
#export FlameManLDFortOpts=""
#export FlameManLDCCOpts=""
#export FlameManRANLIB=""


### change the location of the FlameMaster package here
setenv FlameManPath ~/FlameMaster

### do not change the following couple of lines
setenv FlameManSource $FlameManPath/FlameManSource
setenv FlameManLib $FlameManPath/FlameManLibs
setenv myData $FlameManPath/FlameManData
alias FlameMaster $FlameManLib/FlameMaster
alias LT '$FlameManLib/ListTool -M'
set path=($path $FlameManLib)

### You have to choose between two possible options for FLEXPATH depending
### on weather flex is installed on your machine and in the current path or not.
### You can find out by typing 'which flex' in the command line. If you get a valid
### path as a result, this means flex is installed. Then choose the first option.
### If flex cannot be found, choose the second one. Then, flex will be installed
### during the FlameMaster installation. You will do this following the procedure
### outlined in the Readme file.
###
### Option 1: flex already exists on your machine (don't insert the path in the FLEXPATH definition)
setenv FLEXPATH ""
### Option 2: flex will be installed to $FlameManLib
#setenv FLEXPATH $FlameManLib/bin/

### choose one of the following sets

setenv FlameManMach LINUXGXX
setenv FlameManFF g77
setenv FlameManF90 g90
setenv FlameManCC gcc
setenv FlameManCCOpts "-ansi -O3"
setenv FlameManCXX g++
setenv FlameManCXXOpts "-x c++ -O3"
setenv FlameManF77Opts "-O3"
setenv FlameManLIBF77 ""
setenv FlameManLDFortOpts "-lg2c"
setenv FlameManLDCCOpts ""
setenv FlameManRANLIB 'ranlib $(LIBRARY)'

# following for Mac OSX
#setenv FlameManMach DIGITAL
#setenv FlameManFF g77
#setenv FlameManCC cc
#setenv FlameManCCOpts "-O3"
#setenv FlameManCXX CC
#setenv FlameManCXXOpts "-O3"
#setenv FlameManF77Opts "-O3"
#setenv FlameManLIBF77 "/sw/lib/libg2c.a"
#setenv FlameManLDFortOpts ""
#setenv FlameManLDCCOpts "-lm -lstdc++"
#setenv FlameManRANLIB 'ranlib $(LIBRARY)'

#setenv FlameManMach DIGITAL
#setenv FlameManFF f77
#setenv FlameManCC cc
#setenv FlameManCCOpts ""
#setenv FlameManCXX CC
#setenv FlameManCXXOpts ""
#setenv FlameManF77Opts ""
#setenv FlameManLIBF77 ""
#setenv FlameManLDFortOpts ""
#setenv FlameManLDCCOpts "-lm -lUfor -lfor -lFutil -lots"
#setenv FlameManRANLIB ""

#setenv FlameManMach SUNGXX
#setenv FlameManFF f77
#setenv FlameManCC gcc
#setenv FlameManCCOpts "-ansi"
#setenv FlameManCXX g++
#setenv FlameManCXXOpts "-x c++"
#setenv FlameManF77Opts ""
#setenv FlameManLIBF77 "/optx/SUNWspro/SC4.0/lib/libsunmath.a /optx/SUNWspro/SC4.0/lib/libF77.a /optx/SUNWspro/SC4.0/lib/libm.a"
#setenv FlameManLDFortOpts ""
#setenv FlameManLDCCOpts ""
#setenv FlameManRANLIB 'ranlib $(LIBRARY)'

#setenv FlameManMach SUN
#setenv FlameManFF f77
#setenv FlameManCC cc
#setenv FlameManCCOpts "-Xc"
#setenv FlameManCXX CC
#setenv FlameManCXXOpts ""
#setenv FlameManF77Opts ""
#setenv FlameManLIBF77 "/usr/lang/SC1.0/libF77.a /usr/lang/SC1.0/libm.a"
#setenv FlameManLDFortOpts ""
#setenv FlameManLDCCOpts ""
#setenv FlameManRANLIB ""

# choose following for SGI
#setenv FlameManMach SUN
#setenv FlameManFF f77
#setenv FlameManCC cc
#setenv FlameManCCOpts "-n32"
#setenv FlameManCXX CC
#setenv FlameManCXXOpts "-n32"
#setenv FlameManF77Opts "-n32"
#setenv FlameManLIBF77 ""
#setenv FlameManLDFortOpts "-lftn"
#setenv FlameManLDCCOpts "-n32"
#setenv FlameManRANLIB ""

#setenv FlameManMach HP
#setenv FlameManFF f77
#setenv FlameManCC cc
#setenv FlameManCCOpts "-Aa -c -DHP -z"
#setenv FlameManCXX CC
#setenv FlameManCXXOpts "-Aa -DHP -z +a1"
#setenv FlameManF77Opts "-U"
#setenv FlameManLIBF77 ""
#setenv FlameManLDFortOpts ""
#setenv FlameManLDCCOpts ""
#setenv FlameManRANLIB ""


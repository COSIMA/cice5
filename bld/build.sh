#! /bin/csh -f

set echo on

if ( $#argv < 3 ) then
  echo '*** Please issue the command like ***'
  echo '> ./comp_auscom_cice.sh <platform> <driver> <resolution> [<debug>]'
  echo 'e.g. comp_auscom_cice.sh nci access-om 1440x1080'
  echo 'platform: the machine to run on.'
  echo 'driver: which driver to use.'
  echo 'resolution: grid resolution longitude by latitude.'
  echo 'debug: if this is a unit testing or debug build. Valid options are \'debug\' or \'unit_testing\''
  exit
else
  set platform = $1
  set driver = $2
  set resolution = $3
  set debug = $4
endif

# Location of this model
setenv SRCDIR $cwd
setenv CBLD   $SRCDIR/bld

if ($debug == 'debug') then
    setenv DEBUG 1
else
    setenv DEBUG 0
endif
if ($debug == 'unit_testing') then
    setenv DEBUG 1
    setenv UNIT_TESTING yes
endif

source $CBLD/config.$platform.$driver.$resolution

### Specialty code
setenv CAM_ICE  no        # set to yes for CAM runs (single column)
setenv SHRDIR   csm_share # location of CCSM shared code
setenv IO_TYPE  netcdf    # set to none if netcdf library is unavailable
setenv DITTO    no        # reproducible diagnostics
setenv THRD     no        # set to yes for OpenMP threading
if ( $THRD == 'yes') setenv OMP_NUM_THREADS 2 # positive integer 
setenv BARRIERS yes       # set -Dgather_scatter_barrier, prevents hangs on raijin
setenv AusCOM   yes
if ($driver == 'access') then
    setenv ACCESS   yes
else
    setenv ACCESS   no
endif
setenv OASIS3_MCT yes	  # oasis3-mct version
setenv NICELYR    4       # number of vertical layers in the ice
setenv NSNWLYR    1       # number of vertical layers in the snow
setenv NICECAT    5       # number of ice thickness categories

if ( $AusCOM == 'yes' ) then
    setenv CPLLIBDIR $LIBACCESSOM2_ROOT/build/lib
    setenv CPLLIBS '-L$(CPLLIBDIR)/ -laccessom2'
    setenv CPLINCDIR $LIBACCESSOM2_ROOT/build
    setenv OASISDIR $LIBACCESSOM2_ROOT/oasis3-mct/Linux/build/lib/
    setenv CPL_INCS '-I$(CPLINCDIR)/include -I$(OASISDIR)/psmile.MPI1 -I$(OASISDIR)/mct'
endif

### Setup the version string, this is the git hash of the commit used to build
### the code. The version of an executable can be found with the following
### command: strings <executable> | grep 'CICE_COMMIT_HASH='
setenv GIT_CONFIG_NOGLOBAL 'yes'

set old_hash=`grep 'public :: CICE_COMMIT_HASH =' $SRCDIR/drivers/$driver/version.F90 | cut -d '"' -f 2 | cut -d '=' -f 2`
set new_hash=`git rev-parse HEAD`

if ( $old_hash != $new_hash ) then
    sed -e "s/{CICE_COMMIT_HASH}/$new_hash/g" $SRCDIR/drivers/$driver/version.F90.template > $SRCDIR/drivers/$driver/version.F90
endif

### Location and name of the generated exectuable
setenv EXE cice_${driver}_${resolution}_${NTASK}p.exe

### Where this model is compiled
setenv OBJDIR $SRCDIR/build_${driver}_${resolution}_${NTASK}p
if !(-d $OBJDIR) mkdir -p $OBJDIR

# These variables are set in the appropriate bld/config
@ a = $NXGLOB * $NYGLOB ; @ b = $BLCKX * $BLCKY * $NTASK
@ m = $a / $b ; setenv MXBLCKS $m ; if ($MXBLCKS == 0) setenv MXBLCKS 1
echo Autimatically generated: MXBLCKS = $MXBLCKS

###########################################
# ars599: 24032014
#	copy from /short/p66/ars599/CICE.v5.0/accice.v504_csiro
#	solo_ice_comp
###########################################
### Tracers               # match ice_in tracer_nml to conserve memory
setenv TRAGE   1          # set to 1 for ice age tracer
setenv TRFY    1          # set to 1 for first-year ice area tracer
setenv TRLVL   1          # set to 1 for level and deformed ice tracers
setenv TRPND   1          # set to 1 for melt pond tracers
setenv NTRAERO 0          # number of aerosol tracers 
                          # (up to max_aero in ice_domain_size.F90) 
                          # CESM uses 3 aerosol tracers
setenv TRBRI   0          # set to 1 for brine height tracer
setenv NBGCLYR 7          # number of zbgc layers
setenv TRBGCS  0          # number of skeletal layer bgc tracers 
                          # TRBGCS=0 or 2<=TRBGCS<=9)

### File unit numbers
setenv NUMIN 11           # minimum file unit number
setenv NUMAX 99           # maximum file unit number

if ($IO_TYPE == 'netcdf') then
  setenv IODIR io_netcdf
else if ($IO_TYPE == 'pio') then
  setenv IODIR io_pio
else
  setenv IODIR io_binary
endif



cp -f $CBLD/Makefile.std $CBLD/Makefile

if ($NTASK == 1) then
   setenv COMMDIR serial
else
   setenv COMMDIR mpi
endif
echo COMMDIR: $COMMDIR

setenv DRVDIR $driver

cd $OBJDIR

### List of source code directories (in order of importance).
cat >! Filepath << EOF
$SRCDIR/drivers/$DRVDIR
$SRCDIR/source
$SRCDIR/$COMMDIR
$SRCDIR/$IODIR
$SRCDIR/$SHRDIR
EOF

cc -o makdep $CBLD/makdep.c || exit 2

make VPFILE=Filepath EXEC=$EXE \
           NXGLOB=$NXGLOB NYGLOB=$NYGLOB \
           BLCKX=$BLCKX BLCKY=$BLCKY MXBLCKS=$MXBLCKS \
      -j 8 -f  $CBLD/Makefile MACFILE=$CBLD/Macros.$platform || exit 2

cd ..
pwd
echo NTASK = $NTASK
echo "global N, block_size"
echo "x    $NXGLOB,    $BLCKX"
echo "y    $NYGLOB,    $BLCKY"
echo max_blocks = $MXBLCKS
echo $TRAGE   = TRAGE,   iage tracer
echo $TRFY    = TRFY,    first-year ice tracer
echo $TRLVL   = TRLVL,   level-ice tracers
echo $TRPND   = TRPND,   melt pond tracers
echo $NTRAERO = NTRAERO, number of aerosol tracers
echo $TRBRI   = TRBRI,   brine height tracer
echo $NBGCLYR = NBGCLYR, number of bio grid layers
echo $TRBGCS  = TRBGCS,  number of BGC tracers

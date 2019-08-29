/*
  Debug: Copy write statements from Fortran Code.
         M-x query-replace-regexp in Emacs:
           write(\*,'(a,2f20.14)') 'CV\(.*\) ', \(.*\)
             with
           printf("CV\1 %20.14f %20.14f %20.14f\\n", \2);
           %   with   .
           (\([^[:space:]]*\)))   with   [\1])
           (\([^[:space:]]*\))   with   [\1]
           (\([^[:space:]]*\)\),\([^[:space:]]*\)))   with   [\1][\2])
           (\([^[:space:]]*\),\([^[:space:]]*\))   with   [\1][\2]
           ncl   with   jtot
           ntl   with   jtot3
           nsoil   with   soil.n_soil
           nwiso-1   with   soil.nwater
           tsrf   with   tsfc
           time%    with   time_var
         Delete extra %20.14f or replace with %d for integer and %ld for long.
 */
// Libraries
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

// some differences between gcc and Visual C++ (defines WIN32)
#ifdef WIN32
#include <iostream.h>
#include <fstream.h>
#include <io.h>
#include <direct.h>
#else
#include <sys/stat.h>
#define __max fmax
#define __min fmin
#endif

// Constants
#define code_file "canisotope_v3.4.cpp"            // this file
#define parameter_file "canisotope_v3.4_param.txt" // parameter file
#define sze 42             // canopy layers is 40 with flex
#define sze3 122           // 3 times canopy layers is 120 with flex
#define szeang 19          // number of sky angle classes is 18 with flex
#define soilsze 12         // number of soil layers is 10 with flex
#define c13sze 22          // number of days to remember mean 13C discrimination (+2, flex)
#define watersze 7         // maximum number of isotopic waters + 1 (+2, 0 not used, i.e. flex)
#define betasze 6          // number levels for beta distribution of e.g. lai (+1, 0 not used)
//MC20190702 #define PI 3.1415926536    // pi
#define PI 3.14159265358979323846
//MC20190702 #define PI180 0.0174532925 // pi/180, radians per degree
//MC20190702 #define PI9 2.86478898     // 
//MC20190702 #define PI2 6.2831853072   // 2*pi
#define TN0 273.15         // Celcius <-> Kelvin
#define Gravity 9.81       // Gravity constant
#define vonKarman 0.41     // von Karman's constant
#define waterdensity 1000. // density of water [kg m-3]
#define undef 9999.        // value for output if variable undefined
#define isotolerance 1e-15 // tolerance in isotope variables so that two values are the same

/*
  ----------------------------------------------------

  10-18-2007, Canisotope, revision 3.4

  DENNIS BALDOCCHI
  Ecosystem Science Division
  Department of Environmental Science, Policy and Management
  & Berkeley Atmospheric Science Center
  151 Hilgard Hall
  University of California, Berkeley
  Berkeley, CA 94720-3110
  baldocchi@nature.berkeley.edu
  510-642-2874

  Canisotope is a member of the CANOAK family of models.

  The original isotope version of CANOAK was developed in collaboration with
  Dave Bowling, Universtiy of Utah.
  
  This version includes extensions and updates of Alexander Knohl (ETH Zurich)
  and Matthias Cuntz (MPI-BGC Jena).

  -------------------------------------------------

  CANOAK:

  CANOAK is a coupled biophysical and ecophysiolgical model that
  computes fluxes of water, heat and CO2 exchange within vegetation
  canopies and between the canopy and the atmosphere. In doing so CANOAK
  computes the canopy microclimate (light, wind, temperature, humidity and
  CO2), which provides drivers for physiological processes such as photosynthesis,
  respiration, transpiration and stomatal conductance. The canopy is divided into
  40 layers and at each layer the energy balance, photosynthesis, transpiration,
  stomatal conductance and respiration of sunlit and shaded leaves is computed.
  Stomatal conductance is computed using a model that is linked to photosynthesis.


  The mechanistic and biochemical photosynthesis model of Farquhar is used to scale
  leaf CO2 exchange to the canopy.

  Stomatal conductance is computed using the Ball-Berry routine, that
  is dependent on leaf photosynthesis, relative humidity and the
  CO2 concentration at the leaf's surface.

  Photosynthesis and stomatal conductance are computed using a
  cubic analytical solution I developed from the coupled models
  during my stay in Viterbo, 92/93 (Baldocchi, 1994).

  Photosynthetic parameters (Vcmax and Jmax) are scaled
  with height according to specific leaf weight. SLW is a
  surrogate for the effect of leaf nitrogen on these parameters.
  Kinetic photosynthetic coefficients are derived from Peter
  Harley's 1992 field measurements at WBW. For the seasonal runs Vcmax and
  Jmax are scaled to seasonal changes in leaf area index.

  This model calculates isoprene emissions from the forest canopy using
  algorithms by Guenther.

  Turbulent diffusion is based on Lagrangian theory. A dispersion matrix is
  calculated from the random walk algorithm of Thomson (from program MOVTHM.C).
  The dispersion matrix is scaled according to u*. The dispersion matrix also
  scales with stability, since it is dependent on sigma w and sigma w depends on z/L.
  Dispersion matrices's functional dependence on stability have been derived by
  regression against 5 stability classes.

  Light profiles are computed on the basis of layers with constant
  DELTA Z to be compatible with the meteorological model. Clumping of
  foliage is consider via use of the Markov model for the probability
  of beam penetration. The penetration and scattering routines of Norman (1979) are
  used, these use a slab, one dimensional 'adding' method for scattering.
  The probability of sunlit leaf area is computed using the
  Markov model which introduces a clumping factor to the
  Poisson eq.


  Feedbacks are considered between the source sink function,
  which is a function of C, and the dispersed concentration field,
  which is a function of sources and sinks.

  Solar elevation is computed with the algorithm of Welgraven

  Soil energy balance and heat fluxes are computed using a
  numerical solution to the Fourier heat transfer equation.
  after Campbell.

  Analytical (quadratic) solutions are used for the leaf and soil surface energy balances

  Overall the model requires simple environmental inputs. Solar
  radiation, wind speed, air temperature, humidity, CO2 and a
  deep 32 cm soil temperature, at least, are needed to calculate
  fluxes. Of course parameters describing leaf area, height and
  physiological capacity of the vegetation are needed, but these
  are canopy specific.


  This version computes fluxes over the course of an annual cycle using hourly
  data from Walker Branch Watershed. Inputs include meteorological data and seasonal
  adjustment of leaf area index and photosynthetic capacity. Respiration and electron
  transport rates are scaled with Vcmax, which is a function of leaf nitrogen and depth
  in the canopy.


  Information on the photosynthetic and stomatal conductance model are reported in:

  Baldocchi, D.D. 1994. An analytical solution for coupled leaf photosynthesis
  and stomatal conductance models. Tree Physiology 14: 1069-1079.


  Information on the leaf photosynthetic parameters can be found in:

  Harley, P.C. and Baldocchi, 1995.Scaling carbon dioxide and water vapor exchange
  from leaf to canopy in a deciduous forest:leaf level parameterization.
  Plant, Cell and Environment. 18: 1146-1156.

  Wilson, K.B., D.D. Baldocchi and P.J. Hanson. 2000. Spatial and seasonal variability of
  photosynthesis parameters and their relationship to leaf nitrogen in a deciduous forest.
  Tree Physiology. 20, 565-587.


  Tests of the model are reported in:

  Baldocchi, D.D. 1997. Measuring and modeling carbon dioxide and water vapor
  exchange over a temperate broad-leaved forest during the 1995 summer drought.
  Plant, Cell and Environment. 20: 1108-1122

  Baldocchi, D.D. and P.C. Harley. 1995. Scaling carbon dioxide and water vapor
  exchange from leaf to canopy in a deciduous forest: model testing and application.
  Plant, Cell and Environment. 18: 1157-1173.

  Baldocchi, D.D and T.P. Meyers. 1998. On using eco-physiological, micrometeorological
  and biogeochemical theory to evaluate carbon dioxide, water vapor and gaseous deposition
  fluxes over vegetation. Agricultural and Forest Meteorology 90: 1-26.

  Baldocchi, D.D. Fuentes, J.D., Bowling, D.R, Turnipseed, A.A. Monson, R.K. 1999. Scaling
  isoprene fluxes from leaves to canopies: test cases over a boreal aspen and a mixed species temperate
  forest. J. Applied Meteorology. 38, 885-898.

  Baldocchi, D.D. and K.B.Wilson. 2001. Modeling CO2 and water vapor exchange of a
  temperate broadleaved forest across hourly to decadal time scales. Ecological Modeling 142,155-184

  Baldocchi, D.D. and Bowling, D.R. 2003. Modeling the discrimination of 13CO2 above and within a
  temperate broadleaved forest canopy on hour to seasonal time scales
  Plant, Cell and Environment (in press).

  ----------------------------------------------------

  DEBUGGING NOTES

  October 2007, ak, mc
    - included water cycle
      % interception water
      % soil water (optional)
      % litter water (optional)
    - included water isotopes (optional)
      % ISOLSM - but changed percolation in soil, and changed surface fluxes
    - update 13C-formulation
    - LAI input file (optional)
    - photosynthetic parameter in two canopy layers (optional)
    - new parameter file format
    - only C (no C++), C99 standard; tested on Microsoft Visual C++ and GCC compilers
      the gcc call is: gcc -O3 -lm -x c -std=c99 -o caniso canisotope_v3.4.cpp

  September 2006,ak
  leaf photosynthesis and stomata conductance is now calculated from an analytical solution of the Farquhar model for
  photosynthesis including mesophyll conductance and the Leuning Ball Berry equation for stomata conductance

  August 2006, ak
  diffuse radiation can either be calculated or can come from measurements (input file). Have to check diffuse radiation
  routine since modelled diffuse radiation does not fit well with measurements.

  August 2006, ak
  ustar and z/L are now calculated based on actual measurement height, canopy height, displacement height and roughness length, not fixed parameters

  Jan 2006 ,ak
  soil temperature are now initiatlised with soil temperature from last time step (memory effect)

  Sept 2005, ak
  new routines for soil water added (based on LSM)

  7-6-2004, ak
  Parameter are read in from external file (parameter_file) in same directory. These parameters overwrite default values.
  source code is changed from C to C++ (*.c to *.cpp) to allow object oriented programming.

  12-5-2002
  Changed Dispersion matrix to account for the exponential decrease in sigma w and to use
  the Tl model of Massman and Weil, which considers canopy structure better. It produces a profile
  of Tl that increase with depth, like my data, and produces a relatively constant length scale.
  This makes better sense as L in the canopy is large and scales with the canopy shear length scale,
  eg h, and is representative of large scale coherent eddies.

  7-15-2002. Need to consider diffusion of 13C across the boundary layer, explicitly.

  Found I am using constant diffusivities, irrespective of temperature and pressure, based on
  data from Monteith and Unsworth at 20 C.

  Modified code and adopt reference values of Massman (1998) and am correcting diffusivities for
  temperature and pressure

  added diffuse_par/par as an output

  4-18-2002. Ready to submit paper to journal, but noticed some inconsistancies with R= C13/(C12+C13).
  Also found an error in pdb_CO2. It is 0.01115, not 0.011115. New version has been tested and with
  proper changes, there was only a small change in values. Keeling plot intercept even improved.

  3-28-2002 calling measured CO2 at 44 m on D196
  evaluating del 13C with different functions for day and night
  found algebraic error in recycling calculation

  3-26-2002. I think the signs may be wrong for Fiso. To be consistent with the up and down direction
  of fluxes I am introducing (-)Acan into Fiso

  3-25-2002. Need to average daily Ps for periods when it is not zero. Also computing daily
  ave D13.

  2-28-2002 Found new calcs of recycle are inconsistant with old files. Need to average the layer
  by layer calculations. Now having problems trying to find best way to express it and interpret
  variables in the Yakir and Sternberg Eq.

  2-26- 2002 Sensitivity run with sigma0 =0.3 rather than 0.15 on Dispersion matrix

  2-18-2002 Running with CO2 = 370 ppm based on day values of Bowling, rather than 360

  12-12-2001 Turn off soil respiration; *.resp

  12-6-2001 Sensitivity Run with Vcmax set at 50 umol m-2 s-1

  11/19/2001 Making sure changes to Cansens_v2 have been made to this code.

  found error in Dij function. It should be f=(a+x)/(b+c x)

  6-25-2001 Added the Daaman and Simmonds convection effects to soil resistances.
  Also keeping rh at soil surface too high. Using RH of air to deduce if the soil
  surface may be wet. Then vpdfact=1, else 0.1, like the Canpond simulations.

  Found unit error in the vpd at the soil. Plus it was using the ea at the top of the canopy
  rather than near the soil, as a reference value. Also we need to adjust the vpdsoil factor
  for wet and dry soils. I am assuming the soil is wet if RH is above some level, eg 95%.
  Then vpdfact is one. Otherwise the vpdfact is 0.1, like the test for the Pipo in Oregon.

  6-23-2001. Need to incorporate Daaman and Simmons routines for soil resistances, like
  Canpond. At present soil evaporation is too large during the winter. Thought I had put these
  routines in, but did not.

  6-25-2001. Sorted out the PS root problem and found some errors on the soil evaporation
  routine. Had wrong units and wrong variable for soil Energy balance. Also added Daamen
  and Simmons routine to soil energy. Need to fix this code yet.

  6/17/2001. Fixed minor LAI profiling error
  Incorporated changes in the Photosynthesis routine catch the right roots

  6/7/2001: Be careful. Getting obscenely large Ps values with new root correction

  4/12/2001 Fixed conditional statement in PS().
  Was: if (ROOT1=ROOT2).
  Fixed to: if (ROOT1==ROOT2)


  add recycling to profiles; these need debugging or need to be dropped as
  negative values occur and I am not considering the transfer of air from above
  correctly

  4/11/2001

  Adding computations of recycling of Co2 according to Sternberg

  3/29/2001

  Code is compute photosynthesis weighted D13 and the Keeling int are agreeing with data
  from Bowling.

  Computing fraction of air from soil

  3/23/2001

  Compute Flux weighted D13

  Computing what the assimilated discrimination, R, was to use to respire the next day
  to assess the isotope discrimination of the respired CO2

  3/21/2001

  The current version gives a linear Keeling plot and the intercept is -26 per mill
  the soil value!!


  3/20/2001 Bole.resp was in mg m-2 s-1 and was being added to the source which
  I had converted to umol m-2 s-1. I was also not converting the bole respiration
  to 13C fluxes. This may be one reason Isotopes were getting screwed up in the canopy

  3/13/2001 Redoing the Dij using the algorithm of Massman and Weil. Recomputed
  Dij function for z/L

  Debugging ISOTOPE and CONC

  found that I needed to adjust the CO2 profile by the factor ma/mc rhoair_mole_m-3
  to convert from mole -3 to umol

  Need to make sure 13C and 12C or (12C + 13C) cancel.

  3/12/2001 Got Isotope model running


  2-23-2001 add info on wj and wc profiles

  1-7-2001

  found met.zl was set to zero in CONC for neutral run of dispersion matrix. need to
  turn it off here and in canoak!

  1-2-2001
  Be careful on reading data with C++ compiler. I was getting warning
  errors using floats, so I changed everything to double. When I did this
  the input structure and Dispersion matrices were read in as garbage.
  I had to change those variables back to float and then they read ok!

  NEED TO BE CAREFUL ON ISOTOPES AS MANY OF THE CO2 SOURCE SINKS WERE
  IN UNITS OF MG M-2 S-1. WANT TO BE IN UNITS OF MICROMOLES M-2 S-1

  NEED TO ASSIGN RS TO F(PAR) ONLY THE FIRST ITERATION ON SUN AND
  SHADE LEAVES. AT PRESENT IT IS BEING RESET TO F(PAR) EVERY
  LOOP

  12-28-2000

  converting global variables over to structures

  ----------------------------------------------------
*/

// ----------------------------------------------------
// Declare subroutines and functions
// ----------------------------------------------------

// Radiation transfer routines

void RNET();                     // computes net radiation profile
void PAR();                      // computes visible light profiles
void DIFFUSE_DIRECT_RADIATION(); // computes direct and diffuse radiation from radiation inputs
void NIR();                      // computes near infrared radiation profiles
double SKY_IR(double a);         // computes the incoming longwave radiation flux density
void IRFLUX();                   // computes longwave radiation flux density
void G_FUNC_DIFFUSE();           // computes the leaf orientation function for diffuse radiation
void GFUNC();                    // computes the leaf orientation angle for the given sun angle and a spherical leaf inclination

// canopy structure routines

void ANGLE();            // computes sun angles
double GAMMAF(double a); // gamma function
void LAI_TIME();         // updates LAI with time
void FREQ(double a);     // leaf angle frequency distribution for G

// computes leaf photosynthesis with Farquhar model

void PHOTOSYNTHESIS(double Iphoton,double *rstompt, double zzz, double cca, double tlk,
                    double leleaf, double *A_mgpt, double *GPPpt, double *resppt, double *cipnt, double *cspnt,
                    double *ccpnt, double *cicapnt, double *cccapnt,
                    double *rh_leafpnt, double *vpd_leafpnt, double *wjpnt, double *wcpnt, int JJ);
double TEMP_FUNC(double a,double b,double c,double d, double e); // Arrenhius function for temperature
double TBOLTZ(double a,double b, double c, double d);            // Boltzmann function for temperature
double INV_TBOLTZ(double a, double b, double c, double d);       // inverse of Boltmann function

// computes respirations

void SOIL_RESPIRATION(); // computes soil respiration
void BOLE_RESPIRATION(); // computes bole respiration

// Turbulence and leaf boundary layer routines

double UZ(double a);                                                    // wind speed computation as a function of z
void BOUNDARY_RESISTANCE(double zzz, double TLF, double cws, int JLAY); // leaf boundary layer resistances for heat, water, CO2
void FRICTION_VELOCITY();                                               // updates friction velocity with new z/L

// Concentration calculation routines for q, T and CO2

void CONC(double *ap1, double *ap2, double c, double d, double e);

// leaf energy balance subroutines and functions

void ENERGY_AND_CARBON_FLUXES();       // profiles of energy, stomatal conductance and photosynthesis
void ENERGY_BALANCE(double qrad, double *tsfckpt, double taa, double rhovva, double rvsfc,
                    double stomsfc, double *lept, double *lewet, double *H_leafpt, double *lout_leafpt,
                    double wet_coef); // computes leaf energy balance

double SFC_VPD(double tlk, double Z, double leleafpt); // humidity at leaf surface
double LAMBDA(double a);                               // latent heat of vaporization
double ES(double a);                                   // saturation vapor pressure function
double DESDT(double a);                                // first derivative of saturation vapor pressure function
double DES2DT(double a);                               // second derivative of saturation vapor pressure function

// precipitation throughfall in canopy

void THROUGHFALL(); // calculates throughfall and interception for all layers

// soil energy balance functions and subroutines

void SET_LEAF_PHENOLOGY();       // initialize leaf start and leaf full
void SET_SOIL_TEXTURE();         // set soil texture per soil layer
void SET_SOIL_MOISTURE();        // set humidity per soil layer
void SET_SOIL_ROOT();            // set rooting profile
void SET_SOIL_TEMP();            // initialize deep soil temperature
void SET_SOIL_SAXTON();          // set hydraulic soil parameters after Saxton & Rawls (2006)
void SET_SOIL_CLAPP();           // set hydraulic soil parameters after Clapp & Hornberger (1978)
void SET_SOIL_TIME();            // initialize soil heat conduction equation time step
void SET_SOIL_LITTER_CAPACITY(); // compute soil heat capacity and conductivity
void SOIL_ENERGY_BALANCE();      // soil energy balance
double SOIL_RESISTANCE();
double SOIL_RH();

// declare subroutine for isoprene

void ISOPRENE_LEAF_FLUX(double a, double b, double c, double *dpt); // isoprene emission from leaves
void ISOPRENE_CANOPY_FLUX(double *apt);

// declare subroutine for computing 13C isotope fluxes

void CARBON_ISOTOPE(); // 13C isotope computations

// declare subroutine for input/output

void SKIP_INPUT();
void INPUT_DATA();
void WRITE_DUMMY(int length, ... );

// declare routine for reading in parameter from parameter file

void READ_PARAMETER();

// declare routine for setting parameter and overwriting default values

void SET_PARAMETER();

// house keeping routines

void FILENAMES();      // defines file names and opens files
void INIT_NEXT_STEP(); // zeros arrays & copy former time step

// declare routine for copying files from one location to another

int FileCopy(const char *src, const char *dst );

// declare routine to calculate soil moisture profil

int SOIL_H2O();

// soil moisture helper routines
int TRIDIA();
int IUP(int ihere, double flx);

// declare special water isotope routines

double ISORAT(double *qi, double *q, double *lostqi, double *lostq, int species);
int LE_WISO();
int LEAF_WISO();
int SOIL_FLUX_WISO();
int CANOPY_FLUX_WISO();

// declare litter routines

void SET_LITTER_TEXTURE();
void SET_LITTER_TEMP();
void SET_LITTER_MOISTURE(); // set humidity per soil layer
double LITTER_RESISTANCE();
double LITTER_RH();
void LITTER_RAIN();
void LITTER_H2O();

// declare isotope helper routines

double ALPHA_EQU_H2O(double temp, int species);
double ALPHA_KIN_H2O(double rs, double rb, int species, int merlivat);
double DELTA(double rare, double abundant, double standard);
double DELTA1000(double rare, double abundant, double standard);
double DELTA_H2O(double rare, double abundant, int species);
double DELTA1000_H2O(double rare, double abundant, int species);
double INVDELTA1000_H2O(double iso, int species);
double INVDELTA1000(double iso);
double VTF(double T, int species);

// declare string helper routines

int iswhite(char c);
int find_last_nonwhite_character(char *str);
char *ltrim(char *s);

// declare math helper routines

double corr(int n, double x[], double y[]);

// ----------------------------------------------------
// Declare Global Structures
// ----------------------------------------------------

struct output_variables
{
  double c1;
  double c2;
  double c3;
  double c4;
  double c5;
  double c6;
  double c7;
  double c8;
  double c9;
  double c10;
  double netrad;
  double sumrn;
  double sumlout;
  double sumh;
  double sumle;
  double can_ps_mol;
  double can_gpp;
  double c_transpiration_mole;
  double canresp;
  double sumksi;
  double tleaf_mean;
  double tavg_sun;
  double tavg_shd;
  double ave_D13;
  double ave_D13_long;
  double ave_cica;
  double ave_gs;
  double ave_daC13;
  double diff_par;
  double rnet_soil;
  double dLAIdz[sze];
};
struct output_variables output;

struct flux_variables
{
  double photosyn; // photosynthesis, um m-2 s-1
  double nee; // net ecosystem exchange, umol m-2 s-1
  double c_evaporation[watersze]; // canopy evaporation from interception reservoir, mm m-2 s-1
  double c_transpiration[watersze]; // canopy transpiration through stomata, mm m-2 s-1
  double c_evapotranspiration[watersze]; // canopy transpiration + evaporation, mm m-2 s-1
  double soilevap[watersze]; // soil evaporation, mm m-2 s-1
  double litterevap[watersze]; // litter evaporation, mm m-2 s-1
  double s_evap[watersze]; // soil+litter evaporation, mm m-2 s-1
  double evapotranspiration[watersze]; // canopy+soil+litter evaporationtranspiration, mm m-2 s-1
  double soilinfl[watersze]; // soil infiltration, mm m-2 s-1
  double surfrun[watersze]; // surface runoff = water flux lateral to surface, mm m-2 s-1
  double albedo;				// surface albedo = (PARout + NIRout)/(PARin + NIRin)
};
struct flux_variables flux;

//MC20190615 use doubel precision for input
/* struct input_variables */
/* { */
/*   int dayy; // day */
/*   float hhrr; // hour */
/*   float ta; // air temperature, C */
/*   float rglobal; // global radiation, W m-2 */
/*   float parin; // photosynthetically active radiation, micromole m-2 s-1 */
/*   float pardif; // diffuse PAR, micromol m-2 s-1 */
/*   float lai; // for Nate McDowell's juniper site read lai instead of pardif */
/*   float lai_up; // for Nate McDowell's juniper site read lai in upper canopy */
/*   float lai_down; // for Nate McDowell's juniper site read lai in understory */
/*   float ea; // vapor pressure, kPa */
/*   float wnd; // wind speed, m s-1 */
/*   double ppt[watersze]; // precipitation, mm per hour */
/*   double dppt[watersze]; // delta value of precipitation [per mi] */
/*   double dvapour[watersze]; // delta value of vapour [per mi] */
/*   float co2air; // CO2 concentration, ppm */
/*   float press_mb; // air pressure, mb */
/*   float tsoil; // soil temperature in 50 cm, degree C */
/*   float soilmoisture; // soil moisture in 15 cm, % */
/*   long int flag; // input coding */
/*   float longwave;		// long wave irradiannce, Wm-2 */
/*   float d13CO2;  // d13C of atmospheric CO2 [permille] */
/*   float d18CO2;  // d18O of atmospheric CO2 [permille] */
/* }; */
struct input_variables
{
  int dayy; // day
  double hhrr; // hour
  double ta; // air temperature, C
  double rglobal; // global radiation, W m-2
  double parin; // photosynthetically active radiation, micromole m-2 s-1
  double pardif; // diffuse PAR, micromol m-2 s-1
  double lai; // for Nate McDowell's juniper site read lai instead of pardif
  double lai_up; // for Nate McDowell's juniper site read lai in upper canopy
  double lai_down; // for Nate McDowell's juniper site read lai in understory
  double ea; // vapor pressure, kPa
  double wnd; // wind speed, m s-1
  double ppt[watersze]; // precipitation, mm per hour
  double dppt[watersze]; // delta value of precipitation [per mi]
  double dvapour[watersze]; // delta value of vapour [per mi]
  double co2air; // CO2 concentration, ppm
  double press_mb; // air pressure, mb
  double tsoil; // soil temperature in 50 cm, degree C
  double soilmoisture; // soil moisture in 15 cm, %
  long int flag; // input coding
  double longwave;		// long wave irradiannce, Wm-2
  double d13CO2;  // d13C of atmospheric CO2 [permille]
  double d18CO2;  // d18O of atmospheric CO2 [permille]
};
struct input_variables input;


// structure for time variables
struct time_variables
{
  double local_time;
  long int daytime; // day+hour coding, e.g. 1230830
  int doy;  //day of year
  int days; // day
  int jdold; // previous day
  int count; // number of iterations

  int year0; //always first year of data
  int year; // actual year
  int leafout; // day of leaf out
  int leaffull; // date of full leaf
  int leaffall; // date of leaf shedding
  int leaffallcomplete; // date of leaf shedding completed
  double lai; // lai as a function of time
  double time_step; // time between subsequent computation times
};
struct time_variables time_var;


// structure for meteorological variables
struct meteorology
{
  double ustar; // friction velocity, m s-1
  double ustar_filter; // updated friction velocity with new H, m s-1
  double rhova_g; // absolute humidity, g m-3
  double rhova_kg; // absolute humidity, kg m-3
  double H; // sensible heat flux, W M-2
  double H_filter; // old sensible heat flux, W m-2
  double air_density; // air density, kg m-3
  double T_Kelvin; // absolute air temperature, K
  double zl; // z over L, height normalized by the Obukhov length
  double press_kpa; // station pressure, kPa
  double press_bars; // station pressure, bars
  double press_Pa; // pressure, Pa
  double pstat273; // gas constant computations
  double air_density_mole; // air density, mole m-3
  double relative_humidity; // relative humidity
  
  double dispersion[sze3][sze]; // Lagrangian dispersion matrix, s m-1
};
struct meteorology met;


// structure for surface resistances and conductances
struct surface_resistances
{
  double gcut; // cuticle conductances, mol m-2 s-1
  double fdrought; // drought response factor of stomata conducatance [-]
  double rcuticle[sze]; // cuticle resistance, m2 s mol-1

  double fthreshold; // threshold of relative soil water below which stomata show response to water stress [mm]
  double fslope; // slope of water stress response
};
struct surface_resistances sfc_res;


// structure for plant and physical factors
struct factors
{
  double latent; // latent heat of vaporization, J kg-1
  double heatcoef; // factor for sensible heat flux density
  double a_filt; // filter coefficients
  double co2; // CO2 factor, ma/mc * rhoa (mole m-3)
};
struct factors fact;


// structure for bole respiration and structure
struct bole_respiration_structure
{
  double factor; // base respiration, micromoles m-2 s-1, data from Edwards
  double respiration_mole; // bole respiration, micromol m-2 s-1
  double respiration_mg; // bole respiration, mg CO2 m-2 s-1
  double calc; // calculation factor
  double layer[sze]; // bole pai per layer
};
struct bole_respiration_structure bole;


// structure for canopy architecture
struct canopy_architecture
{
  double bdens[10]; // probability density of leaf angle
};
struct canopy_architecture canopy;


// structure for non dimensional variables
struct non_dimensional_variables
{
  // Prandtl Number

  double pr;
  double pr33;


  // Schmidt number for vapor

  double sc;
  double sc33;

  // Schmidt number for CO2

  double scc;
  double scc33;


  // Schmidt number for ozone

  double sco3;
  double sco333;

  // Grasshof number

  double grasshof;


  // multiplication factors with leaf length and diffusivity

  double lfddv;
  double lfddh;
};
struct non_dimensional_variables non_dim;


// boundary layer resistances
struct boundary_layer_resistances
{
  double vapor; // resistance for water vapor, s/m
  double heat; // resistance for heat, s/m
  double co2; // resistance for CO2, s/m
};
struct boundary_layer_resistances bound_layer_res;


// radiation variables, visible, near infrared and infrared
struct solar_radiation_variables
{
  // profiles of the probabilities of sun and shade leaves
  double prob_beam[sze]; // probability of beam or sunlit fraction
  double prob_shd[sze]; // probability of shade

  double ir_dn[sze]; // downward directed infrared radiation, W m-2
  double ir_up[sze]; // upward directed infrared radiation. W m-2


  // inputs of near infrared radiation and components, W m-2

  double nir_beam; // incoming beam component near infrared radiation, W m-2
  double nir_diffuse; // incoming diffuse component near infrared radiation, W m-2
  double nir_total; // incoming total near infrared radiaion, W m-2

  // computed profiles of near infrared radiation, W m-2

  double nir_dn[sze]; // downward scattered near infrared radiation
  double nir_up[sze]; // upward scattered near infrared radiation
  double nir_sun[sze]; // near infrared radiation on sunlit fraction of layer
  double nir_shd[sze]; // near infrared radiation on shaded fraction of layer
  double beam_flux_nir[sze]; // flux density of direct near infrared radiation

  // inputs of visible light, PAR, W m-2

  double par_diffuse; // diffuse component of incoming PAR, parin
  double par_beam; // beam component of incoming PAR, parin

  // computed profiles of visible radiation, PAR, W m-2

  double par_shd[sze]; // PAR on shaded fraction of layer
  double par_sun[sze]; // PAR on sunlit fraction of layer, beam and diffuse
  double beam_flux_par[sze]; // PAR in the beam component
  double par_down[sze]; // upward scattered PAR
  double par_up[sze]; // downward scattered PAR

  // flux densities of visible quanta on sun and shade leaves for photosynthesis
  // calculations, micromoles m-2 s-1

  double quantum_sun[sze]; // par on sunlit leaves
  double quantum_shd[sze]; // par on shaded leaves

  // Net radiation profiles, W m-2

  double rnet_sun[sze]; // net radiation flux density on sunlit fraction of layer
  double rnet_shd[sze]; // net radiation flux density on shade fraction of layer

  double beta_rad; // solar elevation angle, radians
  double sine_beta; // sine of beta
  double beta_deg; // solar elevation angle, degrees
  double ratrad; // radiation ratio to detect cloud amount

  // leaf and soil optical properities of near infrared radiation

  double nir_reflect; // leaf reflectance in the near infrared
  double nir_trans; // leaf transmittance in the near infrared
  double nir_soil_refl_dry; // soil reflectance in the near infrared of dry soil
  double nir_soil_refl; // soil reflectance in the near infrared
  double nir_absorbed; // leaf absorptance in the near infrared

  // optical properties of leaves and soil for PAR

  double par_absorbed; // PAR leaf absorptance
  double par_reflect; // PAR leaf reflectance
  double par_trans; // PAR leaf transmittance
  double par_soil_refl_dry; // PAR soil reflectance of dry soil
  double par_soil_refl; // PAR soil reflectance

  // Net radiation profiles, W m-2

  double exxpdir[sze]; // exponential transmittance of diffuse radiation through a layer
  double ratradnoon; // radiation ratio at noon for guestimating cloud amount at night
};
struct solar_radiation_variables solar;


// physical properties of soil abd soil energy balance variables
struct soil_variables
{
  // soil properties
  double T_soil[soilsze]; // soil temperature [deg C]
  double T_soil_filter[soilsze]; // soil temperature [deg C]
  double k_conductivity_soil[soilsze]; // thermal conductivity of soil *dz
  double cp_soil[soilsze]; // specific heat of soil, f(texture, moisture)
  double T_Kelvin; // soil surface temperature in Kelvin
  double T_air; // air temperature above soil, C
  double tsfc; // soil surface temperature in C
  double tsfc_filter; // filtered soil surface temperature in C
  double Temp_ref; // reference soil temperature, annual mean, C
  double Temp_amp; // amplitude of soil temperature, C
  double rs; // soil resistances
  double rb; // soil boundary layer resistances
  double T_15cm; // soil temperature at 15 cm

  // soil energy flux densities, W m-2
  double lout; // longwave efflux from soil
  double soilevap; // soil evaporation
  double soilevap_filter; // filtered soil evaporation
  double litterevap; // litter evaporation
  double litterevap_filter; // filtered litter evaporation
  double litterevap_save; // save litter evaporation for wiso
  double maxlitterevap; // maximum possible litter evaporation
  double evap; // total forest floor evaporation (soil+litter)
  double heat; // soil sensible heat flux density
  double rnet; // net radiation budget of the soil
  double gsoil; // soil heat flux density

  // soil CO2 respiratory efflux
  double respiration_mole; // soil respiration, micromol m-2 s-1
  double respiration_mg; // soil respiration, mg CO2 m-2 s-1
  double base_respiration; // base rate of soil respiration, micromol m-2 s-1
  double respiration_auto; // autotrohic soil respiration, micromol m-2 s-1
  double respiration_hetero; // heterotrophic soil respiration, micromol m-2 s-1
  double resp_13; // respiration of 13C micromole m-2 s-1
  double SR_ref_temp; // refernce temperature for soil respiration, degC

  // from ISOLSM
  double swp[soilsze]; //soil water potential (kPa)
  double swp_mm[soilsze]; //soil water potential (mm)
  double r[soilsze][watersze];
  double a[soilsze];
  double b[soilsze];
  double c[soilsze];
  double dwat[soilsze][watersze];

  double theta[soilsze][watersze]; //volumetric soil water content
  double qdrai[watersze]; //sub-surface runoff (mm h2o s-1)
  double qinfl[watersze]; //infiltration rate, mm h2o s-1
  double soil_mm; //water content in the soil, mm m-2
  double soil_mm_root; //total water content in the soil weighted by roots, mm
  double soil_mm_50; //water content in the top 50 cm of soil, mm m-2
  double qtran[watersze]; //plant transpiration, mm s-1
  double qseva[watersze]; //plant transpiration, mm s-1

  // Saxton & Rawls (2006)
  double psi[soilsze]; // soil water potential [Pa]
  double k_theta[soilsze]; // moisture conductivity [mm s-1]

  // Litter (Ogee & Brunet 2002)
  double theta_l[watersze]; // volumetric litter water content
  double T_l; // litter temperature [K]
  double T_l_filter; // litter temperature [K]
  double rl; // litter resistances
  double c_litterevap; // factor regulating that litter evaporation < litter moisture
  double c_litterevap_save; // save factor regulating that litter evaporation < litter moisture
  int maxlitter; // switch if litterevap was restricted to maximum litter evaporation

  double latent; // latent heat coefficient before soil-air vapour pressure deficit
  double lecoef; // latent heat coefficient before soil-air vapour pressure deficit
  double rh_soil; // relative humidity in first soil layer air space
  double rh_litter; // relative humidity in litter layer air space
  
  double xylem[watersze]; // xylem water isotopic composition
  double rxylem[watersze]; // xylem water isotope ratio

  // Fixed variables
  // soil properties
  int n_soil; // # of calculated soil layer
  double z_soil[soilsze]; // depth of soil layer boundaries
  double bulk_density[soilsze]; // soil bulk density
  double d; // dispacement height for soil [m]
  double z0; // soil roughness length [m]
  double T_base; // base soil temperature

  // For numerical soil temperature
  double temperature_dt; // time step, s
  long int temperature_mtime; // number of time steps per dt

  // from ISOLSM
  double moisture_dt; // time step, s
  long int moisture_mtime; // number of time steps per dt
  double qinfl_max; // maximal infiltration rate into soil, mm m-2 s-1
  double clay_in[soilsze]; //clay fraction in each layer input
  double sand_in[soilsze]; //sand fraction in each layer input
  double root[soilsze]; // actual (optimised) relative root abundance (0 to 1)
  double clay[soilsze]; // actual (optimised) clay fraction in each layer
  double sand[soilsze]; // actual (optimised) sand fraction in each layer
  double om[soilsze]; // soil organic matter content in each layer
  double gravel[soilsze]; // volumetric gravel content in each layer
  double root_factor; // gamma of Jackson et al. (Oecologia 1996)
  double sand_factor; // gamma of Jackson et al. (Oecologia 1996)
  double clay_factor; // gamma of Jackson et al. (Oecologia 1996)
  double soil_mm_33_root; //total water content in the soil at field capacity (-33 kPa) weighted by roots, mm
  double soil_mm_1500_root; //total water content in the soil at wilting point (-1500 kPa) weighted by roots, mm

  // Saxton & Rawls (2006)
  int saxton; // 0: Clapp & Hornberger (1978); 1: Saxton & Rawls (2006)
  double theta_1500[soilsze]; // volumetric soil water content at wilting point
  double theta_33[soilsze]; // volumetric soil water content at field capacity
  double theta_s33[soilsze]; // saturation - field capacity volumetric soil water content
  double theta_s[soilsze]; // saturated volumetric soil water content
  double psi_e[soilsze]; // soil water potential at air entry (bubbling pressure) [Pa]
  double rho[soilsze]; // normal soil density (without gravel) [kg m-3]
  double k_s[soilsze]; // saturated conductivity [mm s-1]
  double lambda[soilsze]; // coefficient for moisture conductivity
  double big_a[soilsze]; // coefficient for moisture conductivity
  double big_b[soilsze]; // coefficient for moisture conductivity

  // Litter (Ogee & Brunet 2002)
  double z_litter; // litter layer depth [m]
  double theta_ls; // saturated volumetric litter water content
  double theta_l33; // volumetric litter water content at field capacity
  double n_l; // exponent in litter resistance calculation
  double rho_l; // litter bulk density [kg m-3]

  // Original Dennis code
  int camillo; // 0: Passerat (1986); 1: Camillo & Guerney, i.e. no soil water

  // Water
  double watmin[soilsze]; // residual vol. soil water content depending on texture

  // Drainage
  int drain0; // Switch to write out drainage <0 exactly once
  // Lose water
  int lost0; // Switch to write out when losing water starts
};
struct soil_variables soil;


// Structure for Profile information,fluxes and concentrations
struct profile
{
  // microclimate profiles

  double tair[sze3]; // air temp (C)
  double tair_filter[sze3]; // numerical filter of Tair
  double tair_filter_save[sze3]; // save for later calculations
  double u[sze3]; // wind speed (m/s)
  double rhov_air[sze3][watersze]; // water vapor density for all water isotopes [kg m-3]
  double rhov_air_save[sze3]; // save water vapor density [kg m-3]
  double rhov_air_filter[sze3][watersze]; // numerical filter of rhov_air
  double rhov_air_filter_save[sze3]; // save numerical filter of rhov_air
  double co2_air[sze3]; // co2 concentration (ppm)
  double co2_air_filter[sze3]; // co2 concentration (ppm)
  double vpd_air[sze3]; // vapor pressure deficit[kPa]

  // canopy structure profiles

  double Gfunc_solar[sze];             // leaf-sun direction cosine function
  double tleaf[sze];                   // leaf temperature per layer (sun and shade)
  double isopreneflux[sze];            // isoprene flux per layer (sun and shade)
  double vcmax[sze];                   // vcmax in per layer (temperature corrected)
  double jmax[sze];                    // jmax per layer (temperature corrected)
  double rd[sze];                      // rd per layer
  double vcmaxz[sze];                  // vcmaxz per layer
  double throughfall[sze+1][watersze]; // throughfall in layer [mm], one extra layer for input as precipitation
  double cws[sze][watersze];           // canopy water storage in layer [mm]
  double wet_coef[sze];                // coefficient to regulate evaporation from wet leaf surface [0-1]
  double wet_coef_filter[sze];         // filter coefficient to regulate evaporation from wet leaf surface [0-1]

  // variables for 13C isotopes

  double c13cnc[sze3];       // concentration of 13C
  double c12cnc[sze3];       // concentration of 12C
  double sour13co2[sze];    // source/sink strength 13C
  double d13C[sze3];         // del 13C
  double D13C[sze];          // photosynthetic weighted discrimination, 13D
  double D13C_long[sze];     // photosynthetic weighted discrimination, 13D
  double D13C_ab[sze];	   // boundary layer discrimination
  double D13C_a[sze];		   // stomata discrimination
  double D13C_asal[sze];	   // mesophyll discrimination
  double D13C_b[sze];	       // rubisco discrimination
  double cs[sze];			   // CO2 concentration at leaf surface
  double ci[sze];			   // CO2 concentration inside stomatal cavity
  double cc[sze];			   // CO2 concentration at site of carboxlylation
  double csca[sze];          // ratio of cs and ca
  double cica[sze];          // ratio of ci and ca
  double ccca[sze];          // ratio of cc and ca
  double gm[sze];	     // mesophyll conductance
  double d13Cair[sze3];      // del 13C of the air
  double d13Cplant[sze];     // del 13C of the plant
  double R13_12_air[sze3];   // 13C/12C ratio of the air
  double Rplant_sun[sze];    // ratio of discriminated 13C in sunlit leaves
  double Rplant_shd[sze];    // ratio of discriminated 13C in shaded leaves
  double Rresp[sze];         // discriminated 13C for later respiration, Ps weight
  double Rresp_sum[sze];     // summing of Ps weighted values for daily ave
  double Rresp_ave[sze];     // previous days discriminated 13C for later respiration
  int cnt_Rresp[sze];
  double recycle[sze3];      // fraction of recycled CO2, after Yakir and Sternberg

  // source/sink strengths

  double source_co2[sze]; // source/sink strength of CO2

  // fluxes for total layer, with sun/shade fractions

  double dPsdz[sze];        // layer photosynthesis (micromol m-2 s-1)
  double dPsdz_sun[sze];    // layer photosynthesis of sunlit area (micromol m-2 s-1) = Aphoto sun
  double dPsdz_shd[sze];    // layer photosynthesis of shaded area (micromol m-2 s-1) = Aphoto shd
  double dGPPdz[sze];	      // layer gross primary productivity (Ps + Resp) (micromol m-2 s-1)
  double dGPPdz_sun[sze];	  // layer gross primary productivity (Ps + Resp) of sunlit area (micromol m-2 s-1)
  double dGPPdz_shd[sze];   // layer gross primary productivity (Ps + Resp) of shaded area (micromol m-2 s-1)
  double dHdz[sze];         // layer sensible heat flux (W m-2)
  double dLEdz[sze][watersze]; // layer latent heat flux
  double dLEdz_sun[sze];    // layer latent heat flux of sunlit area (W m-2)
  double dLEdz_shd[sze];    // layer latent heat flux of shaded area (W m-2)
  double dRNdz[sze];        // layer net radiation flux (W m-2)
  double dLoutdz[sze]; // layer radiation flux
  double dRESPdz[sze];      // layer respiration (micromol m-2 s-1)
  double dRESPdz_sun[sze];  // layer respiration of sunlit area (micromol m-2 s-1) = rd sun
  double dRESPdz_shd[sze];  // layer respiration of shaded area (micromol m-2 s-1) = rd shd
  double dStomCondz[sze];   // layer stomatal conductance (mol m-2 s-1)
  double dStomCondz_sun[sze];   // layer stomatal conductance (mol m-2 s-1)
  double dStomCondz_shd[sze];   // layer stomatal conductance (mol m-2 s-1)

  // sun leaf variables

  double sun_frac[sze]; // sun leaf fraction
  double sun_tleaf[sze]; // leaf temp (C)
  double sun_GPP[sze]; // layer GPP flux for sun only (micromol mn-2 s-1)
  double sun_A[sze]; // layer A flux for sun only (micromol mn-2 s-1)
  double sun_gs[sze]; // stomatal conductance (m s-1)
  double sun_gs_mol[sze]; // stomatal conductance of sun leaves mol m-2 s-1
  double sun_rs[sze]; // stomatal resistance to H2O (s/m)
  double sun_rs_filter[sze]; // filter stomatal resistance to H2O (s/m)
  double sun_rs_save[sze]; // stomatal resistance to H2O (s/m)
  double sun_rbh[sze]; // boundary layer resistance to heat (s/m)
  double sun_rbv[sze]; // boundary layer resistance to H2O (s/m)
  double sun_rbco2[sze]; // boundary layer resistance to CO2 (s/m)
  double sun_cs[sze]; // Cs on sun leaves (CO2 mixing ratio on leaf surface
  double sun_ci[sze]; // Ci on sun leaves
  double sun_cc[sze]; // Cc (CO2 mixing ratio at site of carboxylation) on sun leaves
  double sun_D13[sze]; // discrimination 13C
  double sun_D13_ab[sze]; // del 13C boundary resistence on shaded leaves
  double sun_D13_a[sze]; // del 13C stomatal resistance on shaded leaves
  double sun_D13_asal[sze];// del 13C mesophyll resistance on shaded leaves
  double sun_D13_b[sze]; // del 13C rubisco on shaded leaves
  double sun_D13_long[sze];// discrimination 13C
  double sun_csca[sze]; // Cs/Ca on sunlit leaves
  double sun_cica[sze]; // Ci/Ca on sunlit leaves
  double sun_ccca[sze]; // Cc/Ca on sunlit leaves
  double sun_lai[sze]; // sunlit lai of layer
  double sun_tleaf_filter[sze]; // filtered sunlit temperature
  double sun_wj[sze]; // electron transport rate of Ps for sun leaves
  double sun_wc[sze]; // carboxylatio velocity for sun leaves
  double sun_resp[sze]; // respiration
  double sun_isopreneflux[sze]; // isoprene flux per layer for sunleaves
  double iso_sun[sze]; // isoprene flux per leaf area in the sun
  double sun_vpd[sze]; // vapor pressure deficit at leaf surface [hPa]
  double sun_rh[sze]; // relative humidity at leaf surface [hPa]
  double sun_LEstoma[sze][watersze];            // LE flux through stomata [W m-2 leaf]
  double sun_LEwet[sze][watersze];         // LE flux from wet leaves [W m-2 leaf]

  // shade leaf variables

  double shd_frac[sze]; // shade leaf fraction
  double shd_tleaf[sze]; // temperature of shaded leaves
  double shd_GPP[sze]; // photosynthesis (GPP) of shaded leaves
  double shd_A[sze]; // photosynthesis of shaded leaves
  double shd_gs[sze]; // stomatal conductance of shade leaves
  double shd_gs_mol[sze]; // stomatal conductance of shade leaves mol m-2 s-1
  double shd_rs[sze]; // stomatal resistance of shaded leaves
  double shd_rs_filter[sze]; // filter stomatal resistance of shaded leaves
  double shd_rs_save[sze]; // stomatal resistance of shaded leaves
  double shd_rbh[sze]; // boundary layer resistance for heat on shade leaves
  double shd_rbv[sze]; // boundary layer resistance for vapor on shade leaves
  double shd_rbco2[sze]; // boundary layer resistance for CO2 on shade leaves
  double shd_cs[sze]; // Cs on shade leaves (CO2 mixing ratio on leaf surface
  double shd_ci[sze]; // Ci on shade leaves
  double shd_cc[sze]; // Cc (CO2 mixing ratio at site of carboxylation) on sun leaves
  double shd_D13[sze]; // del 13C on shaded leaves
  double shd_D13_ab[sze]; // del 13C boundary resistence on shaded leaves
  double shd_D13_a[sze]; // del 13C stomatal resistance on shaded leaves
  double shd_D13_asal[sze];// del 13C mesophyll resistance on shaded leaves
  double shd_D13_b[sze]; // del 13C rubisco on shaded leaves
  double shd_D13_long[sze]; // discrimination 13C
  double shd_csca[sze]; // Cs/Ca ratio on shaded leaves
  double shd_cica[sze]; // Ci/Ca ratio on shaded leaves
  double shd_ccca[sze]; // Cc/Ca ratio on shaded leaves
  double shd_lai[sze]; // shaded lai of layer
  double shd_tleaf_filter[sze]; // previous temperature
  double shd_wc[sze]; // carboxylation rate for shade leaves
  double shd_wj[sze]; // electron transport rate for shade leaves
  double shd_resp[sze]; // respiration
  double shd_isopreneflux[sze]; // isoprene flux per layer for shade leaves
  double iso_shd[sze]; // isoprene flux per leaf area in the shade
  double shd_vpd[sze]; // vapor pressure deficit at leaf surface [hPa]
  double shd_rh[sze]; // relative humidity at leaf surface [hPa]
  double shd_LEstoma[sze][watersze]; // shd leaf LE flux through stomata [W m-2 leaf]
  double shd_LEwet[sze][watersze];  // LE flux from wet shaded leaves [W m-2 leaf]

  // water isotope variables

  double sun_wi[sze]; // 
  double shd_wi[sze]; // 
  double wa[sze]; // 
  double sun_h[sze]; // 
  double shd_h[sze]; // 
  double rs_fact[sze]; // 
  double sun_gross[sze]; // 
  double shd_gross[sze]; // 
  double sun_LEstoma_new[sze]; // sun leaf LE flux through stomata in accordance with water isotope [mol(H2O)/m2(leaf)s]
  double shd_LEstoma_new[sze]; // shd leaf LE flux through stomata in accordance with water isotope [mol(H2O)/m2(leaf)s]
  double sun_LEstoma_save[sze]; // save LE flux through stomata [W m-2 leaf]
  double shd_LEstoma_save[sze];   // save LE flux from wet leaves [W m-2 leaf]

  double sun_alpha_k[sze][watersze]; // 
  double shd_alpha_k[sze][watersze]; // 
  double sun_alpha_equ[sze][watersze]; // 
  double shd_alpha_equ[sze][watersze]; // 
  double rvapour[sze3][watersze]; // 
  double sun_peclet[sze][watersze]; // 
  double shd_peclet[sze][watersze]; // 
  double sun_fem[sze][watersze]; // 
  double shd_fem[sze][watersze]; // 
  double sun_craig[sze][watersze]; // 
  double shd_craig[sze][watersze]; // 
  double sun_leafwater_e[sze][watersze]; // 
  double shd_leafwater_e[sze][watersze]; // 
  double sun_leafwater_e_old[sze][watersze]; // 
  double shd_leafwater_e_old[sze][watersze]; // 
  double sun_leafwater[sze][watersze]; // 
  double shd_leafwater[sze][watersze]; // 
  double sun_trans_rtrans[sze][watersze]; // 
  double shd_trans_rtrans[sze][watersze]; // 
  double sun_rtrans[sze][watersze]; // 
  double shd_rtrans[sze][watersze]; // 

  // canopy structure profiles
  double ht[sze]; // layer height (m)
  double dLAIdz[sze]; // leaf area index of layer (m2/m2)
  double dPAIdz[sze]; // plant area index of layer
  double Gfunc_sky[sze][szeang]; // leaf-sky sector direction cosine function
};
struct profile prof;


struct switch_variable
{
  int soil_resp_temp; // 1: temp for hetero resp = input soil; 0: modelled soil at 5cm
  int mode; // 1: optimisation mode
  int ball; // 0: Ball-Berry; 1: Leuning
  int isoprene; // 1: calc isoprene
  int d13c; // 1: calc 13CO2
  int wiso; // 1: calc water isotopes
  int bethy_resp; // 0: auto/hetero=50/50; 1: Bethy
  int no_neg_water_flux; // 1: restrict water fluxes >= 0
};
struct switch_variable set_switch;


struct isotope_variable
{
  double delta_soil; // delta13C of heterotrophic soil respiration
  double da_m_day; // slope of daytime regression deltaCa*Ca=m*Ca+b
  double da_b_day; // intercept of daytime regression deltaCa*Ca=m*Ca+b
  double da_m_night; // slope of nighttime regression deltaCa*Ca=m*Ca+b
  double da_b_night; // intercept of nighttime regression deltaCa*Ca=m*Ca+b
  double bigdelta[c13sze]; // daily flux weighted canopy discrimination after simple equation for 20 days
  double bigdelta_long[c13sze]; // daily flux weighted canopy discrimination after extended equation for 20 days
};
struct isotope_variable Cisotope;


struct water_isotope_variable
{
  int nwater; // 0: if no water isotopes; x:# of isotopic waters + normal water
  int merlivat; // 0: Cappa et al.; 1: Merlivat
  int nofracsoil; // 0: normal; 1: no fractionation at soil evaporation
  int nofraclitter; // 0: normal; 1: no fractionation at litter evaporation
  int nofracleaf; // 0: normal; 1: no fractionation at leaf transpiration
  int nofracin; // 0: normal; 1: rain same as initial
  int implicit; // 0: explicit, i.e. out of loop; 1:implicit, i.e. in loop
  double vsmow[watersze];
  double lost[watersze]; //stores water lost by rounding errors during ratio calculation
  double dtheta[soilsze][watersze]; // initial soil delta-values
};
struct water_isotope_variable wiso;

// ----------------------------------------------------
// Declare Default Parameter Values
// ----------------------------------------------------

// Parameter are here set on default values (plus for gloabl definition).
// Later new parameters are read in from parameter_file in read_parameter()

int nparameter; // # of parameters read in read_parameter()
char parameter_name[500][256]; // declare parameter name array that will be read in from parameter_file in read_parameter()
char parameter_value[500][256]; // declare parameter array that will be read in from parameter_file in read_parameter()

// file and directory structure

char disksubdir[256]; // working directory
char filesuffix[256]; // output suffix
char input_file_subdir[256]; // input file location
char input_file_suffix[256]; // input file suffix
char dispfile[256]; // dispersion file

char outputfolder[128]; // save from filenames() to print out in main

// Site location

double latitude  = 51.07; // Hainich
double longitude = 10.45; // 

// Eastern Standard TIME
double zone = -1.; // delay from GMT

// canopy structure variables

double ht = 33.0; // Canopy height [m]
double zm = 43.0; // measurement height [m]
double zd = 22; // displacement height [m]
double z0 = 3.3; // roughness length [m]

double pai = 1.4; // Plant area index [m2 m-2]
double lai = 6.4; // Leaf area index [m2 m-2]

// model can use either vcmax and jmax at optimal temperature or at 25 deg C
// vcmax(25 deg C) and j__max(25 deg C) is later converted to vcopt and jmopt with INV_TBOLTZ()
double vcopt[sze];// = 80; // carboxylation rate at optimal temperature, umol m-2 s-1
double jmopt[sze];// = 170; // electron transport rate at optimal temperature, umol m-2 s-1
double vc25[sze];// = 66.; // carboxylation rate at 25 deg C, umol m-2 s-1
double jm25[sze];// = 127.5; // electron transport rate at 25 deg C, umol m-2 s-1
double rd25 = 0.34; // dark respiration at 25 C, rd25= 0.34 umol m-2 s-1, currently not used since rd25 is calculated as fraction of vcmax

double jm_vc[sze];// = 2; //ratio of jmax to vcmax : varies between 1.5 and 2 see von Caemmerer 2000
double rd_vc[sze];// = 0.015; //ratio of rd to vcmanx : varies between 0.01 to 0.02 see von Caemmerer 2000

double ustar_ref = 1.00; //reference ustar for Dij

int jtot=40; // number of canopy layers
int jtot1=41; // jtot + 1 = above canopy
int jtot3=120; // number of layers in the domain, three times canopy height
int izref=52; //changed: array value of reference ht at 43 m, zm/ht*jtot

//MC20190702 double pi4=12.5663706;
double PI2   = 2.0 * PI;  // 2*pi
double pi4   = 4.0 * PI;  // 4*pi - not used
double PI9   = 9.0 / PI;  // 9/pi
double PI180 = PI / 180.; // pi/180

double delz = 0.825; // changed: height of each layer, ht/jtot
double zh65=0.019697; // changed: 0.65/ht

// Universal gas constant

double rugc = 8.3144; // J mole-1 K-1
double rgc1000 = 8314.4; // gas constant times 1000.

double Rw = 461.5; // gas constant for water vapour [J kg-1 K-1] (was 461.89 in v2.0)

// Consts for Photosynthesis model and kinetic equations.
// for Vcmax and Jmax. Taken from Harley and Baldocchi (1995, PCE)

double hkin = 200000.0; // enthalpy term
double skin = 710.0; // entropy term
double ejm = 48000.0; // activation energy for electron transport
double evc = 55000.0; // activation energy for carboxylation

// Enzyme constants & partial pressure of O2 and CO2
// Michaelis-Menten K values. From survey of literature.

double kc25 = 274.6; // kinetic coef for CO2 at 25 C, microbars = 27.46 Pa
double ko25 = 419.8; // kinetic coef for O2 at 25C, millibars = 41980 Pa
double o2 = 210.0; // oxygen concentration umol mol-1 (should be mmol mol-1 = 212 hPa = 21200 Pa)

// tau is computed on the basis of the Specificity factor (102.33)
// times Kco2/Kh2o (28.38) to convert for value in solution
// to that based in air/
// The old value was 2321.1.

// New value for Quercus robor from Balaguer et al. 1996

double tau25 = 2904.12; // tau coefficient

// Arrhenius constants
// Eact for Michaelis-Menten const. for KC, KO and dark respiration
// These values are from Harley

double ekc = 80500.0; // Activation energy for K of CO2; J mol-1
double eko = 14500.0; // Activation energy for K of O2
double erd = 38000.0; // activation energy for dark respiration, eg Q10=2
double ektau = -29000.0;
double tk_25 = 298.15; // absolute temperature at 25 C
double toptvc = 311.0; // optimum temperature for maximum carboxylation
double toptjm = 311.0; // optimum temperature for maximum electron transport
double eabole=45162; // activation energy for bole respiration for Q10 = 2.02

// Constants for leaf energy balance

double sigma = 5.67e-08; // Stefan-Boltzmann constant W M-2 K-4
double cp = 1005.; // Specific heat of air, J KG-1 K-1
double mass_air = 28.97; // Molecular mass of dry air, g mole-1
double mass_CO2=44.; // molecular mass of CO2, g mole-1
double mass_13CO2=45.; // molecular mass of 13CO2, g mole-1
double dldt = -2370.; // Derivative of the latent heat of vaporization

double ep = 0.98; // emissivity of leaves
double epsoil = 0.98; // Emissivity of soil
double epsigma = 5.5566e-8; // ep*sigma
double epsigma2 = 11.1132e-8; // 2*ep*sigma
double epsigma4 = 22.2264e-8; // 4.0 * ep * sigma
double epsigma6 = 33.3396e-8; // 6.0 * ep * sigma
double epsigma8 = 44.448e-8; // 8.0 * ep * sigma
double epsigma12= 66.6792e-8; // 12.0 * ep * sigma

double betfact = 1.5; // multiplication factor for aerodynamic
// sheltering, based on work by Grace and Wilson

// constants for the polynomial equation for saturation vapor pressure-T function, es=f(t)
double a1en=617.4;
double a2en=42.22;
double a3en=1.675;
double a4en=0.01408;
double a5en=0.0005818;

// Leuning model for stomatal conductance

double g0[sze];
double a1[sze];
double D0[sze];

// Ball-Berry stomatal coefficient for stomatal conductance

double kball[sze];// = 9.5;

// intercept of Ball-Berry model, mol m-2 s-1

double bprime[sze];// = 0.0175; // intercept for H2O

// Minimum stomatal resistance, s m-1.

double rsm = 145.0;
double brs = 60.0; // curvature coeffient for light response

// Ratio of mesophyll conductance to Vcmax

double gm_vc = 0.008; 

// curvature of light response curve

double curvature = 0.8; 

// leaf quantum yield, electrons

double qalpha = 0.22;
double qalpha2 = 0.0484; // qalpha squared, qalpha2 = pow(qalpha, 2);

// leaf clumping factor

double markov = 0.84;

// Leaf dimension. geometric mean of length and width (m)

double lleaf = 0.1; // leaf length, m
double attfac = 2.5; // attenuation factor for wind speed inside canopy

// LAI beta-distribution

double ht_midpt[betasze]; // mid points of beta distribution integration
double lai_freq[betasze]; // fractions of total lai per layer
double pai_freq[betasze]; // fractions of total pai per layer

// Optical properties of leaves

double par_reflect[5];       // PAR leaf reflectivity (winter: bark reflectivity, avg top and bottom)
double par_trans[5];         // PAR leaf transmissivity
double par_soil_refl_dry[5]; // PAR soil reflectivitiy
double nir_reflect[5];       // NIR leaf reflectivity
double nir_trans[5];         // NIR leaf transmissivity
double nir_soil_refl_dry[5]; // NIR soil reflectivitiy

// Precipitation interception

double water_film_thickness = 0.2; // typically 0.1 [mm m-2 LAI] on each side of leaf
double tau_water = 1.0; // typically 0.5 interception efficiency per leaf area

// Diffusivity values for 273 K and 1013 mb (STP) using values from Massman (1998) Atmos Environment
// These values are for diffusion in air. When used these values must be adjusted for
// temperature and pressure

// nu, Molecular viscosity

double nuvisc = 13.27; // mm2 s-1
double nnu = 0.00001327; // m2 s-1

// Diffusivity of CO2

double dc = 13.81; // mm2 s-1
double ddc = 0.00001381; // m2 s-1

// Diffusivity of heat

double dh = 18.69; // mm2 s-1
double ddh = 0.00001869; // m2 s-1

// Diffusivity of water vapor

double dv = 21.78; // mm2 s-1
double ddv = 0.00002178; // m2 s-1

// Diffusivity of ozone

double do3 = 14.44; // mm2 s-1
double ddo3 = 0.00001444; // m2 s-1

// Diffusivity of 13CO2

// Isotope ratio of PeeDee Belimdite standard (PDB)
// redefined as 13C/(12C+13C) for PDB, Tans et al 93

double Rpdb_CO2 = 0.01115;

// Isotope ratio of PeeDee Belimdite standard (PDB)
// defined as 13C/12C, from Farquhar

double Rpdb_12C = 0.01124;

// Start and End of model run in year, format dddhhmm

long start_run = 0000000; //start on doy 000 at 0000 hours
long end_run = 3660000; //end on doy 366 at 0000 hours

long start_profiles = 1400000; //start model profiles on doy 140 at 0000 hours
long end_profiles = 1450000; //end model profiles on doy 145 at 0000 hours

// leaf development

long leaf_out = 120; //start of leaf out
long leaf_full = 160; //day of full leaves
long leaf_fall = 290; //leaf fall
long leaf_fall_complete = 320; //end of leaf fall

// number of leaf sides with stomata
double n_stomata_sides = 1.;

// Declare file pointers for I/O

// fptr1   met in
// fptr4   dispersion matrix
// fptr6   season out
// fptr7   soil out
// fptr8   daily out
// fptr9   profile out
// fptr10  flux profile out
// fptr11  optimise
// fptr12  dummy out
// fptr13  h2o soil out
// fptr14  water iso in
// fptr15  lai in
// fptr16  wiso h2o leaf out
// fptr17  wiso profile out
// fptr18  wiso flux profile out
// fptr19  wiso h2o soil out
// fptr20  13c daily out
// fptr21  13c season out
// fptr22  13c profile air out
// fptr23  13c profile flux out
FILE *fptr1,*fptr6, *fptr7, *fptr4, *fptr8, *fptr9, *fptr10, *fptr11;
FILE *fptr12, *fptr13, *fptr14, *fptr15, *fptr16, *fptr17, *fptr18, *fptr19;
FILE *fptr20, *fptr21, *fptr22; //, *fptr23; 

// reached last time
int lastin;

// # of lines in input file
// long lines_in_file;

// Nate McDowell extra stuff

int extra_nate = 0; // Do extras for Nate McDowell's New Mexico site
// only used if extra_nate==1
int nup = 9; // 1 to nup-1 is understory, nup to jtot is upper canopy

// ----------------------------------------------------
// Main routine
// ----------------------------------------------------

// ----------------------------------------------------
int main()
{
  int ok=1; // main return

  // sets parameter

  // Declare internal variables

  int i_count=0, i_max=0, i=0, j=0, junk=0;

  long int daycnt=0, cntps=0;

  char header[300]; //header of input file

  //MC20190615 float dummy=0.;
  double dummy=0.;

  double rnet_soil=0., netrad=0.;

  double fc_mg=0., fc_mol=0., wue=0., fc_13C=0.;

  double sun_A=0., shd_A=0., sumA=0.;
  double sumfc=0.,sumevap=0.,sumps=0., sumgpp=0., sumsens=0., sumlai=0.; // sums for daily averages and totals
  double sumpar=0.,sumnet=0.,sumta=0.,sumbole=0., sumsoil=0., sumgs=0.;
  double sumh=0., sumle=0., sumrn=0., sumlout=0., sumresp=0., sumksi=0., sum13C=0., sumbolelay=0., sumF13C=0.;
  double ebalance=0., sbalance=0., tbalance=0.;
  double sumD13=0., sumD13_long=0., sumD13_long_isotope=0.;
  double ave_D13=0., ave_D13_long=0., ave_D13_long_isotope=0., ave_daC13=0.;
  double sum_csca=0., sum_cica=0., sum_ccca=0., sum_gs=0., sum_gm=0.;
  double ave_csca=0., ave_cica=0., ave_ccca=0., ave_gs=0., ave_gm=0.;
  double sumD13C=0., sumD13C_long_day=0., sumTleaf=0.;
  double sumD13_a=0., sumD13_ab=0., sumD13_asal=0., sumD13_b=0.;
  double ave_D13_a=0., ave_D13_ab=0., ave_D13_asal=0., ave_D13_b=0.;

  double R_air_12C=0.,  Rlongterm_12C=0;

  double can_ps_mol=0., can_gpp=0., can_ps_mg=0., canresp=0., c_transpiration_mole=0.;

  //MC double da13=0., Fisotope=0., ave_Ca_f=0., ave_d13resp=0.;

  double isoprene_efflux=0.; // , recycle=0.;

  double etest=0., etest_old1=0., etest_old2=0., etest_diff1=0., etest_diff2=0.; // test value for energy closure
  double itest=0., itest_old1=0., itest_old2=0., itest_diff1=0., itest_diff2=0.; // test for wiso convergence

  double tleaf_mean=0., tavg_sun=0., tavg_shd=0.; // integrated leaf temperatures

  int mc=0;
  double tmp1[sze3], tmp2[sze3]; // arrays for data transfer with CONC()
  double tmp3=0., temp3=0.;
  clock_t start, end; // time model runs
  double elapsed=0.;
  int doit=0, count40=0, totalcount=0, dcount40=0, dtotalcount=0, ncount40=0, ntotalcount=0;
  double rcws[watersze], ea=0.;
  double temp1=0.;

  int print_balance=0;

  // To have same header as Fortran for easier debug
  /* printf("\n\n\n\n\n\n\n\n\n\n\n\n\n"); */
  printf("\nStart Canoak.\n\n");

  start = clock();

  // Isotope standards
  
  // Pee Dee Belimdite (PDB) for 13CO2
  Rpdb_12C = 0.01124; // 13C/12C from Farquhar
  Rpdb_CO2 = 0.01115; // redefined as 13C/(12C+13C), Tans et al. (1993)

  // Vienna Standard Mean Ocean Water (VSMOW) for water isotopes
  wiso.vsmow[2] = 2005.2e-6;  // 18O/16O
  wiso.vsmow[3] = 155.76e-6;  // 2H/1H
  wiso.vsmow[4] = 1;          // virtual H2O

  soil.d = 1e-6; // displacement height for soil = 0.
  soil.z0 = 0.015; // roughness length of soil = 1.5 cm
  soil.drain0 = 0;
  soil.lost0 = 0;
  i_max = 101; // maximum number of iterations for energy balance
  
  READ_PARAMETER(); // reads in parameter array from parameter_file
  SET_PARAMETER();

  // soil.T_base is set to mean annual air temperature
  if (extra_nate == 1)
    soil.T_base = 10.8; // Nate McDowells Juniper site
  else
    soil.T_base = 7.5; // Hainich
  
  if (soil.n_soil == 0)
    soil.camillo = 1; // 0: layered soil water; 1: no soil water model

  // Water isotopes or not: nwater is 0 if no water isotopes
  if (set_switch.wiso == 1)
    wiso.nwater = 3;
  else
    wiso.nwater = 0;

  // Input year to name input files

  time_var.year = time_var.year0;

  // 'secure' T profile
  soil.T_base = __min(__max(soil.T_base, -10.), 50.);

  // Define Filenames, Constructs Filenames and Opens files

  FILENAMES();

  SET_LEAF_PHENOLOGY(); // set leaf onset and full
  /* printf("CV01 %d %d\n", time_var.leafout, time_var.leaffull); */
  /* printf("CV02 %d %d\n", time_var.leaffall, time_var.leaffallcomplete); */

  SET_SOIL_TEXTURE(); // set soil texture per soil layer

  SET_SOIL_ROOT(); // set rooting profile

  SET_SOIL_MOISTURE(); // set humidity per soil layer

  SET_SOIL_TEMP(); // set deep soil temperature
  /* printf("CV03 %20.14f %20.14f %20.14f\n", soil.theta_ls, soil.theta_l33, soil.n_l); */
  /* printf("CV04 %20.14f %20.14f\n", soil.root[1], soil.root[soil.n_soil]); */
  /* printf("CV05 %20.14f %20.14f\n", soil.theta[1][1], soil.theta[soil.n_soil][1]); */
  /* printf("CV06 %20.14f %20.14f\n", soil.T_soil[0], soil.T_soil[soil.n_soil+1]); */
  /* printf("CV07 %20.14f %20.14f\n", soil.T_soil_filter[0], soil.T_soil_filter[soil.n_soil+1]); */

  SET_LITTER_TEXTURE();

  SET_LITTER_TEMP();

  SET_LITTER_MOISTURE();
  /* printf("CV08 %20.14f %20.14f %20.14f\n", soil.theta_ls, soil.theta_l33, soil.n_l); */
  /* printf("CV09 %20.14f %20.14f\n", soil.T_l, soil.T_l_filter); */
  /* printf("CV10 %20.14f\n", soil.theta_l[1]); */
  
  if (soil.saxton == 1) // set hydraulic soil parameters
    {
      SET_SOIL_SAXTON();
      /* printf("CV11.01 %20.14f %20.14f\n", soil.theta_1500[1], soil.theta_1500[soil.n_soil]); */
      /* printf("CV11.02 %20.14f %20.14f\n", soil.theta_33[1], soil.theta_33[soil.n_soil]); */
      /* printf("CV11.03 %20.14f %20.14f\n", soil.theta_s33[1], soil.theta_s33[soil.n_soil]); */
      /* printf("CV11.04 %20.14f %20.14f\n", soil.psi_e[1], soil.psi_e[soil.n_soil]); */
      /* printf("CV11.05 %20.14f %20.14f\n", soil.theta_s[1], soil.theta_s[soil.n_soil]); */
      /* printf("CV11.06 %20.14f %20.14f\n", soil.rho[1], soil.rho[soil.n_soil]); */
      /* printf("CV11.07 %20.14f %20.14f\n", soil.big_b[1], soil.big_b[soil.n_soil]); */
      /* printf("CV11.08 %20.14f %20.14f\n", soil.big_a[1], soil.big_a[soil.n_soil]); */
      /* printf("CV11.09 %20.14f %20.14f\n", soil.lambda[1], soil.lambda[soil.n_soil]); */
      /* printf("CV11.10 %20.14f %20.14f\n", soil.k_s[1], soil.k_s[soil.n_soil]); */
      /* printf("CV11.11 %20.14f\n", soil.soil_mm_33_root); */
      /* printf("CV11.12 %20.14f\n", soil.soil_mm_1500_root); */
      /* printf("\n"); */
    }
  else
    {
      SET_SOIL_CLAPP();
      /* printf("CV12.01 %20.14f %20.14f\n", soil.theta_s[1], soil.theta_s[soil.n_soil]); */
      /* printf("CV12.02 %20.14f %20.14f\n", soil.k_s[1], soil.k_s[soil.n_soil]); */
      /* printf("CV12.03 %20.14f %20.14f\n", soil.psi_e[1], soil.psi_e[soil.n_soil]); */
      /* printf("CV12.04 %20.14f %20.14f\n", soil.big_b[1], soil.big_b[soil.n_soil]); */
      /* printf("CV12.05 %20.14f %20.14f\n", soil.theta_1500[1], soil.theta_1500[soil.n_soil]); */
      /* printf("CV12.06 %20.14f %20.14f\n", soil.theta_33[1], soil.theta_33[soil.n_soil]); */
      /* printf("CV12.07 %20.14f %20.14f\n", soil.theta_s33[1], soil.theta_s33[soil.n_soil]); */
      /* printf("CV12.08 %20.14f %20.14f\n", soil.rho[1], soil.rho[soil.n_soil]); */
      /* printf("CV12.09 %20.14f %20.14f\n", soil.k_s[1], soil.k_s[soil.n_soil]); */
      /* printf("CV12.10 %20.14f\n", soil.soil_mm_33_root); */
      /* printf("CV12.11 %20.14f\n", soil.soil_mm_1500_root); */
    }

  // calculate relative plant available water
  soil.soil_mm_root = 0.;
  for (j=1; j<=soil.n_soil; j++)
    soil.soil_mm_root += soil.root[j]*soil.theta[j][1]*1000.
      *(soil.z_soil[j]-soil.z_soil[j-1])*(1.-soil.gravel[j]/100.);
  /* printf("CV13 %20.14f\n", soil.soil_mm_root); */

  // initialize some variables

  // Constants for leaf boundary layers

  non_dim.lfddh = lleaf/ddh;

  // Prandtl Number

  non_dim.pr = nuvisc / dh;
  non_dim.pr33 = pow(non_dim.pr, 0.33);

  // DIFFUSIVITY OF WATER VAPOR, m2 s-1

  non_dim.lfddv = lleaf/ddv;

  // SCHMIDT NUMBER FOR VAPOR

  non_dim.sc = nuvisc / dv;
  non_dim.sc33 = pow(non_dim.sc, 0.33);

  // SCHMIDT NUMBER FOR CO2

  non_dim.scc = nuvisc / dc;
  non_dim.scc33 = pow(non_dim.scc, 0.33);

  // Grasshof Number

  //non_dim.grasshof = Gravity*pow(lleaf, 3)/pow(nnu, 2);
  non_dim.grasshof = Gravity*lleaf*lleaf*lleaf/nnu/nnu;

  // boundary layer resistances

  bound_layer_res.heat = 0.;
  bound_layer_res.vapor = 0.;
  bound_layer_res.co2 = 0.;
  
  for (j=1; j<=jtot3; j++)
    {
      tmp1[j] = 0.;
      tmp2[j] = 0.;
    }

  // assign heights to array

  for (j=1; j<=jtot; j++)
    prof.ht[j]= delz * (double) j;
  /* printf("CV14 %20.14f %20.14f %20.14f\n", non_dim.lfddh, non_dim.pr, non_dim.pr33); */
  /* printf("CV15 %20.14f %20.14f %20.14f\n", non_dim.lfddv, non_dim.sc, non_dim.sc33); */
  /* printf("CV16 %20.14f %20.14f %20.14f\n", non_dim.scc, non_dim.scc33, non_dim.grasshof); */
  /* printf("CV17 %20.14f %20.14f\n", prof.ht[1], prof.ht[jtot]); */

  /* *********************************************************************
     MAIN PROGRAM

     Describe canopy attributes

     Flow chart of the Main Program:

     1a) start run at beginning of year
     1b) establish new LAI
     1c) Compute continuous distrib of height at a function of lai
     1d) Set soil parameters

     2) input met values

     3) Compute solar elevation angle;

     4) Compute PAR and NIR profiles and their flux densities on
     the sunlit and shaded leaf fractions

     5) Compute sunlit and shaded leaf fractions

     6) Compute first estimate of stomatal conductance

     7) Compute first estimate of IRFLUX assuming leaf temperature
     equals air temperature.
     8) Compute first estimate of leaf energy balance and leaf
     temperature

     9) Compute photosynthesis, transpiration

     10) Compute soil energy balance

     11) update computation of friction velocity with new H and z/L

     12) Compute new scalar and source/sink profiles for CO2, T and e

     13) Iterate among 7 through 12 until convergence

     14) compute fluxes of isoprene, 13C isotopes or whatever

     15) Thats all Folks!!!

     ******************************************************************* */

  /* define parameters for soil energy balance model. This version only
     computes heat transfer. A water transfer module needs to be added.
     Preliminary tests show that the soil model can be sensitive to the depth
     of the litter layer
  */

  SET_SOIL_TIME();
  /* printf("CV18 %20.14f %ld\n", soil.temperature_dt, soil.temperature_mtime); */
  /* printf("CV19 %20.14f %20.14f %ld\n", soil.moisture_dt, soil.moisture_dt, soil.moisture_mtime); */

  // Input data on Thomson dispersion matrix that was computed offline with
  // MOVOAK.C, Dij (s m-1)
  for (j=1; j<=jtot; j++)
    {
      for (i=1; i<=jtot3; i++)
        {
          //MC20190615 fscanf(fptr4,"%f %i\n", &dummy, &junk);
          fscanf(fptr4,"%lf %i\n", &dummy, &junk);
          met.dispersion[i][j] = dummy;
        } // next i
    } // next j
  /* printf("CV20 %20.14f %20.14f %20.14f\n", met.dispersion[1][1], met.dispersion[jtot3/2][jtot/2], met.dispersion[jtot3][jtot]); */

  // Initialize respired 13C to the heterotrophic soil respiration value 
  // (assumed to equal longterm average), then update each day

  if (set_switch.d13c == 1)
    {
      Rlongterm_12C = INVDELTA1000(Cisotope.delta_soil)*Rpdb_12C;
      for (j=1; j<=jtot; j++)
	{
	  prof.Rresp_ave[j] = Rlongterm_12C;
	  prof.Rresp_sum[j] = 0.;
	  prof.cnt_Rresp[j] = 0;
	}
  
      for (j=1; j<=20; j++)  //just some arbitary values are taken for initialisation
	{
	  Cisotope.bigdelta[j]      = 18.;
	  Cisotope.bigdelta_long[j] = 18.;
	}
     /* printf("CV21 %20.14f\n", prof.Rresp_ave[1]); */
    }

  // loop through the input met file. There should be a line of data for
  // each hour of the year

  // read and skip first line of input file (header)
  fgets(header, 300, fptr1);
  if (set_switch.wiso == 1) fgets(header, 300, fptr14);
  if (extra_nate == 1) fgets(header, 300, fptr15); // lai file
  SKIP_INPUT(); // skip lines of input file before start_run

  time_var.days  = start_run / 10000;
  time_var.jdold = start_run / 10000;
  // set input.dayy because of several years at once
  input.dayy = time_var.days; // for more then 1 year
  /* printf("CV22 %d %d %d\n", time_var.days, time_var.jdold, input.dayy); */
  if (print_balance == 0)
    {
      printf("Day-of-year\n%03i,", time_var.days);
      if ((time_var.days % 24) == 1) printf("\n");
    }

  // initialise arbitrary value
  solar.ratradnoon = 0.5;

  // initialise interception reservoir
  for (mc=1; mc<=wiso.nwater+1; mc++)
    for (j=1; j<=jtot; j++)
      prof.cws[j][mc] = 0.;

  // initialise leaf water at evaporating site
  for (mc=2; mc<=wiso.nwater+1; mc++)
    for (j=1; j<=jtot; j++)
      {
	prof.sun_leafwater_e[j][mc]     = wiso.vsmow[mc];
	prof.shd_leafwater_e[j][mc]     = wiso.vsmow[mc];
	prof.sun_leafwater_e_old[j][mc] = wiso.vsmow[mc];
	prof.shd_leafwater_e_old[j][mc] = wiso.vsmow[mc];
      }
    
  // initialise leaf area index
  //MC20190615 use better precision for lai
  /* input.lai = (float) lai; */
  /* if (extra_nate == 1) */
  /*   { */
  /*     // for Nate McDowell's juniper site read LAI instead of diffuse PAR */
  /*     input.lai_up   = (float) lai; */
  /*     input.lai_down = (float) lai; */
  /*   } */
  input.lai = lai;
  if (extra_nate == 1)
    {
      // for Nate McDowell's juniper site read LAI instead of diffuse PAR
      input.lai_up   = lai;
      input.lai_down = lai;
    }
  LAI_TIME(); // define leaf area and canopy structure
  /* printf("CV23 %20.14f\n", time_var.lai); */
  /* printf("CV24 %20.14f %20.14f %20.14f\n", solar.par_reflect, solar.par_trans, solar.par_soil_refl_dry); */
  /* printf("CV25 %20.14f %20.14f %20.14f\n", solar.par_absorbed, solar.nir_reflect, solar.nir_trans); */
  /* printf("CV26 %20.14f %20.14f\n", solar.nir_soil_refl_dry, solar.nir_absorbed); */
  /* printf("CV27 %20.14f %20.14f\n", prof.dLAIdz[1], prof.dLAIdz[jtot]); */
  /* printf("CV28 %20.14f %20.14f\n", prof.dPAIdz[1], prof.dPAIdz[jtot]); */
  /* printf("CV29 %20.14f %20.14f\n", solar.exxpdir[1], solar.exxpdir[jtot]); */
  
  // Initialize humidity and temperature profiles, arbitrary values
  for (j=1; j<=jtot; j++)
    {
      prof.sun_tleaf[j]        = 25.;
      prof.shd_tleaf[j]        = 25.;
      prof.sun_tleaf_filter[j] = 25.;
      prof.shd_tleaf_filter[j] = 25.;
      prof.sun_rs[j]           = 1./(bprime[j]*rugc*TN0/101325.);
      prof.sun_rs_filter[j]    = 1./(bprime[j]*rugc*TN0/101325.);
      prof.shd_rs[j]           = 1./(bprime[j]*rugc*TN0/101325.);
      prof.shd_rs_filter[j]    = 1./(bprime[j]*rugc*TN0/101325.);
    }
  /* printf("CV30 %20.14f %20.14f\n", prof.sun_rs[1], prof.sun_rs[jtot]); */

  for (j=1; j<=jtot3; j++)
    {
      prof.tair[j]                 = 25.;
      prof.tair_filter[j]          = 25.;
      prof.tair_filter_save[j]     = 25.;
      prof.rhov_air[j][1]          = 0.02;
      prof.rhov_air_save[j]        = 0.02;
      prof.rhov_air_filter[j][1]   = 0.02;
      prof.rhov_air_filter_save[j] = 0.02;
      prof.co2_air[j]              = 380.;
      prof.co2_air_filter[j]       = 380.;
    }
  if (set_switch.d13c == 1)
    {
      R_air_12C = INVDELTA1000(-8.)*Rpdb_12C; // ratio of 13C relative to 12C
      for (j=1; j<=jtot3; j++)
	{
	  prof.R13_12_air[j] = R_air_12C;
	  prof.d13Cair[j]    = -8.;
	}
    }
  /* printf("CV31 %20.14f %20.14f\n", prof.R13_12_air[1], prof.R13_12_air[jtot3]); */

  for (mc=2; mc<=wiso.nwater+1; mc++)
    for (j=1; j<=jtot3; j++)
      {
	prof.rhov_air[j][mc]        = 0.02*wiso.vsmow[mc];
	prof.rhov_air_filter[j][mc] = 0.02*wiso.vsmow[mc];
	prof.rvapour[j][mc]         = wiso.vsmow[mc];
      }
  /* printf("CV32 %20.14f %20.14f\n", prof.rhov_air[1][2], prof.rhov_air[jtot3][3]); */

  for (j=1; j<=jtot; j++)
    {
      prof.dHdz[j]       = 0.0;
      prof.source_co2[j] = 0.0;
    }
  //MC20190615 for (mc=1; mc<=wiso.nwater; mc++)
  for (mc=1; mc<=wiso.nwater+1; mc++)
    for (j=1; j<=jtot; j++)
      prof.dLEdz[j][mc] = 0.0;
  
  // initialize soil surface temperature with air temperature
  soil.tsfc        = 25.;
  soil.tsfc_filter = 25.;

  // initialise evapotranspiration
  soil.soilevap          = 0.;
  soil.soilevap_filter   = 0.;
  soil.litterevap        = 0.;
  soil.litterevap_filter = 0.;

  // initialise sensible heat flux
  met.H        = 0.;
  met.H_filter = 0.;

  // Start time steps
  doit = 1;
  count40 = 0;
  totalcount = 0;
  dcount40 = 0;
  dtotalcount = 0;
  ncount40 = 0;
  ntotalcount = 0;
  while (doit)
    {
      INPUT_DATA(); // input met data

      INIT_NEXT_STEP(); // zero variables and copy former time steps
      /* printf("CV33.01 %d %20.14f %20.14f\n", input.dayy, input.hhrr, input.ta); */
      /* printf("CV33.02 %20.14f %20.14f %20.14f\n", input.ea, input.wnd, input.ppt[1]); */
      /* printf("CV33.03 %20.14f %20.14f %20.14f\n", input.co2air, input.press_mb, input.tsoil); */
      /* printf("CV33.04 %20.14f %ld %20.14f\n", input.soilmoisture, input.flag, input.d13CO2); */
      /* printf("CV33.05 %20.14f %d %ld\n", input.d18CO2, time_var.year, time_var.daytime); */
      /* printf("CV33.06 %d %20.14f %d\n", time_var.doy, time_var.local_time, time_var.days); */
      /* printf("CV33.07 %20.14f %20.14f %20.14f\n", met.T_Kelvin, met.rhova_g, met.rhova_kg); */
      /* printf("CV33.08 %20.14f %20.14f %20.14f\n", met.press_kpa, met.press_bars, met.press_Pa); */
      /* printf("CV33.09 %20.14f %20.14f\n", met.relative_humidity, met.pstat273); */
      /* printf("CV33.10 %20.14f %20.14f\n", sfc_res.rcuticle[1], sfc_res.rcuticle[jtot]); */
      /* printf("CV33.11 %20.14f %20.14f\n", input.co2air, input.parin); */
      /* printf("CV33.13 %20.14f %20.14f\n", input.wnd, met.ustar_filter); */
      /* printf("CV33.14 %20.14f %20.14f %20.14f\n", input.ta, met.air_density, met.air_density_mole); */
      /* printf("CV33.15 %20.14f %20.14f %20.14f\n", input.dppt[1], input.dppt[2], input.dppt[3]); */
      /* printf("CV33.16 %20.14f %20.14f %20.14f\n", input.dppt[4], input.dvapour[2], input.dppt[3]); */
      /* printf("CV33.17 %20.14f %20.14f\n", wiso.dtheta[1][2], wiso.dtheta[1][wiso.nwater]); */
      /* printf("CV33.18 %20.14f %20.14f\n", input.ppt[2], input.ppt[wiso.nwater]); */
      /* printf("CV33.19 %20.14f %20.14f %20.14f\n", input.lai_up, input.lai_down, input.lai); */
      /* printf("CV33.20 %20.14f %20.14f %20.14f\n", input.rglobal, input.parin, input.pardif); */

      if (time_var.daytime > end_run)
	{
	  doit = 0;
	  goto endday; // short circuit
	}
      if (input.flag == 1) goto endday;

      /* Breakpoints
	 if (time_var.days == 110.) puts("stop");
	 if (time_var.daytime == 1860150) 
	 {
	 puts("stop");
	 }
      */

      // if new day then compute daily averages and update LAI

      // use input.dayy instead of time_var.days because of several years at once
      if (input.dayy > time_var.jdold)
	{
	  // compute daily averages
	  if (cntps == 0)
	    {
	      sumgpp = 0;
	      sumTleaf = undef;
	    }
	  else
	    {
	      sumgpp /= cntps;
	      sumTleaf /= cntps;
	    }
	  if (set_switch.d13c == 1)
	    {
	      if (cntps == 0)
		{
		  sumD13C = undef;
		  sumD13C_long_day = undef;
		}
	      else
		{
		  sumD13C /= sumps;
		  sumD13C_long_day /=sumps;
		}
	    }
	  sumfc /= daycnt;
	  sumevap /= daycnt;
	  sumsens /= daycnt;
	  sumpar /= daycnt;
	  sumnet /= daycnt;
	  sumbole /= daycnt;
	  sumsoil /= daycnt;
	  sumta /= daycnt;
	  sumgs /= daycnt;
	  sumresp /= daycnt;
	  sumps /= daycnt;
	  if (set_switch.d13c == 1)
	    sumF13C /= daycnt;

	  if (set_switch.d13c == 1)
	    {
	      // compute the previous days isotope ratio based on its photosynthesis
	      if (sumps > 0)
		for (j=1; j<=jtot; j++)
		  prof.Rresp_ave[j] = prof.Rresp_sum[j]/prof.cnt_Rresp[j];
	      else
		for (j=1; j<=jtot; j++)
		  prof.Rresp_ave[j] = Rlongterm_12C;

	      for (j=1; j<=jtot; j++)
		{
		  prof.Rresp_sum[j] = 0;
		  prof.cnt_Rresp[j] = 0;
		}

	      // shift memory bigdelta one day back
	      if (sumps > 0)
		{
		  for (j=1; j<=20; j++)
		    {
		      Cisotope.bigdelta[j+1] = Cisotope.bigdelta[j];
		      Cisotope.bigdelta[1] = sumD13C;
		      Cisotope.bigdelta_long[j+1] = Cisotope.bigdelta_long[j];
		      Cisotope.bigdelta_long[1] = sumD13C_long_day;
		    }
		}
	      else  // in case of no photosynthesis keep bigdelta from previous day
		{
		  for (j=1; j<=20; j++)
		    {
		      Cisotope.bigdelta[j+1] = Cisotope.bigdelta[j];
		      Cisotope.bigdelta[1] = Cisotope.bigdelta[j+1];
		      Cisotope.bigdelta_long[j+1] = Cisotope.bigdelta_long[j];
		      Cisotope.bigdelta_long[1] = Cisotope.bigdelta_long[j+1];
		    }
		}  
	    }

	  // Daily output file
	  fprintf(fptr8,"%i,%f,%f,%10.1f,%10.1f,%10.1f,%10.1f,"
		  "%g,%10.2f,%10.2f,%10.2f,%10.2f,%14.10f,%10.3f\n",
		  time_var.jdold, sumfc, sumevap, sumsens, sumpar, sumnet, time_var.lai,
		  sumps, sumresp, sumbole, sumsoil, sumta, sumgs, sumTleaf);
	  if (set_switch.d13c == 1)
	    fprintf(fptr20,"%i,%10.4f,%10.3f\n",
		    time_var.jdold, sumF13C, sumD13C);

	  if (extra_nate == 0) // done below at each time step for Nate McDowells juniper site
	    LAI_TIME(); // update LAI with new day

	  // Re-zero summing variables

	  sumfc = 0.;
	  sumevap = 0.;
	  sumsens = 0.;
	  sumps = 0.;
	  sumgpp = 0.;
	  sumpar = 0.;
	  sumnet = 0.;
	  sumbole = 0.;
	  sumsoil = 0.;
	  sumta = 0.;
	  sumgs = 0.;
	  sumresp = 0.;
	  sumTleaf = 0.;
	  sumF13C = 0.;
	  sumD13C = 0.;
	  sumD13C_long_day = 0.;

	  cntps = 0;
	  daycnt = 0;

	  if (print_balance == 0)
	    {
	      printf("%03i,", time_var.days);
	      if ((time_var.days % 24) == 1) printf("\n");
	    }
	} // end if new day

	  // for Nate McDowell's juniper site read LAI instead of diffuse PAR
	  //   update LAI in each time step
      if (extra_nate == 1)
	LAI_TIME();
      /* printf("CV34.01 %20.14f\n", time_var.lai); */
      /* printf("CV34.02 %20.14f %20.14f %20.14f\n", solar.par_reflect, solar.par_trans, solar.par_soil_refl_dry); */
      /* printf("CV34.03 %20.14f %20.14f %20.14f\n", solar.par_absorbed, solar.nir_reflect, solar.nir_trans); */
      /* printf("CV34.04 %20.14f %20.14f\n", solar.nir_soil_refl_dry, solar.nir_absorbed); */
      /* printf("CV34.05 %20.14f %20.14f\n", prof.dLAIdz[1], prof.dLAIdz[jtot]); */
      /* printf("CV34.06 %20.14f %20.14f\n", prof.dPAIdz[1], prof.dPAIdz[jtot]); */
      /* printf("CV34.07 %20.14f %20.14f\n", solar.exxpdir[1], solar.exxpdir[jtot]); */

      // Compute solar elevation angle

      ANGLE();
      /* printf("CV35 %20.14f %20.14f %20.14f\n", solar.beta_rad, solar.sine_beta, solar.beta_deg); */

      // make sure PAR is zero at night. Some data have negative offset, 
      //   which causes numerical problems

      if (solar.sine_beta <= 0.05)
	input.parin = 0.;
      if (input.parin < 0.)
	input.parin = 0.;

      // Compute the fractions of beam and diffuse radiation from incoming measurements

      // Set the radiation factor to the day before for night calculations. This way if the day
      //   was cloudy, so will IR calculations for the night reflect this.
      if (solar.sine_beta > 0.05)
	{
	  DIFFUSE_DIRECT_RADIATION();
	}
      else
	{
	  solar.ratrad = solar.ratradnoon;
	  //MC20190625 Moved from input_data()
	  solar.par_beam = 0.;
	  solar.par_diffuse = 0.;
	  solar.nir_beam = 0.;
	  solar.nir_diffuse = 0.;
	}
      /* printf("CV36.01 %20.14f %20.14f %20.14f\n", solar.ratrad, output.c10, solar.ratradnoon); */
      /* printf("CV36.02 %20.14f %20.14f\n", solar.par_beam, solar.par_diffuse); */
      /* printf("CV36.03 %20.14f %20.14f\n", solar.nir_beam, solar.nir_diffuse); */

      // computes leaf inclination angle distribution function, the mean direction cosine
      // between the sun zenith angle and the angle normal to the mean leaf

      // for CANOAK we use the leaf inclination angle data of Hutchison et al. 1983, J Ecology
      // for other canopies we assume the leaf angle distribution is spherical

      if (solar.sine_beta > 0.05)
	GFUNC();
      else
	for (j=1; j<=jtot; j++)
	  prof.Gfunc_solar[j] = 0.01;
      /* printf("CV37 %20.14f %20.14f\n", prof.Gfunc_solar[1], prof.Gfunc_solar[jtot]); */

      // Set soil reflectivity depending on soil moisture of first layer
      //   i.e. wet soil seems have as bright as dry soil
      //   after Wilson & Henderson-Sellers (1985)
      //   cited in Knorr PhD (1997, Table 2.4, p. 34)
      tmp3 = __min(__max((soil.theta[1][1]-soil.watmin[1])
			 /(soil.theta_s[1]-soil.watmin[1]), 0.), 1.);
      solar.par_soil_refl = solar.par_soil_refl_dry * (1.-0.5*tmp3);
      solar.nir_soil_refl = solar.nir_soil_refl_dry * (1.-0.5*tmp3);
      /* printf("CV38.01 %20.14f %20.14f\n", solar.par_soil_refl_dry, solar.nir_soil_refl_dry); */
      /* printf("CV38.02 %20.14f %20.14f\n", solar.par_soil_refl, solar.nir_soil_refl); */
      /* printf("CV38.03 %20.14f %20.14f %20.14f\n", soil.theta[1][1], soil.watmin[1], soil.theta_s[1]); */
      
      // Compute PAR profiles

      PAR();
      /* printf("CV39.01 %20.14f %20.14f\n", prof.sun_lai[1], prof.sun_lai[jtot]); */
      /* printf("CV39.02 %20.14f %20.14f\n", prof.shd_lai[1], prof.shd_lai[jtot]); */
      /* printf("CV39.03 %20.14f %20.14f\n", solar.prob_beam[1], solar.prob_beam[jtot]); */
      /* printf("CV39.04 %20.14f %20.14f\n", solar.prob_shd[1], solar.prob_shd[jtot]); */
      /* printf("CV39.05 %20.14f %20.14f\n", solar.par_down[1], solar.par_down[jtot]); */
      /* printf("CV39.06 %20.14f %20.14f\n", solar.par_up[1], solar.par_up[jtot]); */
      /* printf("CV39.07 %20.14f %20.14f\n", solar.beam_flux_par[1], solar.beam_flux_par[jtot]); */
      /* printf("CV39.08 %20.14f %20.14f\n", solar.quantum_shd[1], solar.quantum_shd[jtot]); */
      /* printf("CV39.09 %20.14f %20.14f\n", solar.quantum_sun[1], solar.quantum_sun[jtot]); */
      /* printf("CV39.10 %20.14f %20.14f\n", solar.par_shd[1], solar.par_shd[jtot]); */
      /* printf("CV39.11 %20.14f %20.14f\n", solar.par_sun[1], solar.par_sun[jtot]); */

      // Compute NIR profiles

      NIR();
      /* printf("CV40.01 %20.14f %20.14f\n", solar.nir_up[1], solar.nir_up[jtot]); */
      /* printf("CV40.02 %20.14f %20.14f\n", solar.nir_dn[1], solar.nir_dn[jtot]); */
      /* printf("CV40.03 %20.14f %20.14f\n", solar.nir_shd[1], solar.nir_shd[jtot]); */
      /* printf("CV40.04 %20.14f %20.14f\n", solar.nir_sun[1], solar.nir_sun[jtot]); */
      /* printf("CV40.05 %20.14f %20.14f\n", solar.beam_flux_nir[1], solar.beam_flux_nir[jtot]); */
      /* printf("CV40.06 %20.14f\n", solar.nir_total); */

      // Interception reservoir

      THROUGHFALL();
      /* printf("CV41.01 %20.14f %20.14f\n", prof.throughfall[1][1], prof.throughfall[jtot][wiso.nwater]); */
      /* printf("CV41.02 %20.14f %20.14f\n", prof.cws[1][1], prof.cws[jtot][wiso.nwater]); */
      /* printf("CV41.03 %20.14f %20.14f\n", prof.wet_coef_filter[1], prof.wet_coef_filter[jtot]); */

      // Update litter with rain

      if (soil.camillo == 0 && soil.z_litter >= 1e-6)
	LITTER_RAIN();

      // calculate soil heat capacity and conductivity based on soil moisture, texture and bulk density

      SET_SOIL_LITTER_CAPACITY();

      fact.heatcoef = met.air_density * cp;

      // any recursive looping would occur here, below TL initialization

      i_count = 0;
      time_var.count = 0;
      etest_old1 = 0.;
      etest_old2 = -1.;
      itest_old1 = 0.;
      itest_old2 = -1.;
      fact.a_filt = 0.85;

      met.ustar = met.ustar_filter;

      if (soil.z_litter < 1e-6)
	soil.c_litterevap = 0.;
      else
	soil.c_litterevap = 1.;
      /* printf("CV42.01 %20.14f %20.14f\n", soil.cp_soil[1], soil.cp_soil[soil.n_soil]); */
      /* printf("CV42.02 %20.14f %20.14f\n", soil.k_conductivity_soil[1], soil.k_conductivity_soil[soil.n_soil]); */
      /* printf("CV42.03 %20.14f %20.14f\n", fact.heatcoef, met.ustar); */

      // iteration looping for energy fluxes and scalar fields
      // iterate until energy balance closure occurs or 75 iterations are
      // reached
      do
	{
	  // save input filtered variables for use in wiso routines
	  for (j=1; j<=jtot3; j++)
	    {
	      prof.tair_filter_save[j] = prof.tair_filter[j];
	      prof.rhov_air_filter_save[j] = prof.rhov_air_filter[j][1];
	      prof.rhov_air_save[j] = prof.rhov_air[j][1];
	    }
          /* printf("CV43.01 %20.14f %20.14f\n", prof.tair_filter_save[1], prof.tair_filter_save[jtot3]); */
          /* printf("CV43.01 %20.14f %20.14f\n", prof.rhov_air_filter_save[1], prof.rhov_air_filter_save[jtot3]); */
          /* printf("CV43.01 %20.14f %20.14f\n", prof.rhov_air_save[1], prof.rhov_air_save[jtot3]); */
	  
	  IRFLUX(); // first guess
	  /* printf("CV44.01 %20.14f %20.14f\n", solar.ir_up[1], solar.ir_up[jtot]); */
	  /* printf("CV44.02 %20.14f %20.14f\n", solar.ir_dn[1], solar.ir_dn[jtot]); */

	  FRICTION_VELOCITY();
          /* printf("CV45 %20.14f %20.14f\n", met.zl, met.ustar); */

	  // compute net radiation balance on sunlit and shaded leaves
	  RNET();
	  /* printf("CV46.01 %i %20.14f %20.14f\n", i_count, solar.par_sun[1], solar.nir_sun[1]); */
	  /* printf("CV46.02 %20.14f %20.14f\n", solar.par_sun[40], solar.nir_sun[40]); */
	  /* printf("CV46.03 %20.14f %20.14f\n", solar.par_shd[1], solar.nir_shd[1]); */
	  /* printf("CV46.04 %20.14f %20.14f\n", solar.par_shd[40], solar.nir_shd[40]); */
	  /* printf("CV46.05 %20.14f %20.14f\n", solar.ir_dn[1], solar.ir_up[1]); */
	  /* printf("CV46.06 %20.14f %20.14f\n", solar.ir_dn[40], solar.ir_up[40]); */
	  /* printf("CV46.07 %20.14f %20.14f\n", solar.rnet_sun[1], solar.rnet_shd[1]); */
	  /* printf("CV46.08 %20.14f %20.14f\n", solar.rnet_sun[40], solar.rnet_shd[40]); */

	  // Compute leaf energy balance, leaf temperature, photosynthesis and stomatal conductance.
	  ENERGY_AND_CARBON_FLUXES();
	  /* printf("CV47.01 %i %20.14f %20.14f\n", i_count, prof.dLEdz[1][1], prof.dLEdz[40][1]); */
	  /* printf("CV47.02 %20.14f %20.14f\n", prof.dHdz[1], prof.dHdz[40]); */
	  /* printf("CV47.03 %20.14f %20.14f\n", prof.dRNdz[1], prof.dRNdz[40]); */
	  /* printf("CV47.04 %20.14f %20.14f\n", prof.dLoutdz[1], prof.dLoutdz[40]); */
	  /* printf("CV47.00 %20.14f %20.14f\n", prof.sun_tleaf[1], prof.sun_tleaf[40]); */
	  /* printf("CV47.05 %20.14f %20.14f\n", prof.shd_tleaf[1], prof.shd_tleaf[40]); */

	  // Soil energy balance
	  SOIL_ENERGY_BALANCE();
	  /* printf("CV47.06 %20.14f %20.14f %20.14f\n", flux.soilevap[1], flux.litterevap[1], flux.s_evap[1]); */
	  /* printf("CV47.07 %20.14f %20.14f %20.14f\n", soil.tsfc, soil.T_l, soil.litterevap); */
	  /* printf("CV47.08 %20.14f %20.14f %20.14f\n", output.c1, output.c2, output.c3); */
	  /* printf("CV47.09 %20.14f %20.14f\n", output.c4, output.c7); */
	  /* printf("CV47.10 %20.14f %20.14f %20.14f\n", soil.theta[1][1], soil.theta[2][1], soil.theta[3][1]); */

	  // --- Water ---

	  if (set_switch.wiso == 1)
	    {
	      // calculates latent heat in accordance with leaf water isotope calculations
	      LE_WISO();

	      // re-define/re-calc LEstoma in accordance with leaf water isotope calculations
	      for (j=1; j<=jtot; j++)
		{
		  prof.sun_LEstoma_save[j] = prof.sun_LEstoma[j][1];
		  prof.shd_LEstoma_save[j] = prof.shd_LEstoma[j][1];
		  prof.sun_LEstoma[j][1] = prof.sun_LEstoma_new[j]*LAMBDA(prof.tair_filter[j]+TN0);
		  prof.shd_LEstoma[j][1] = prof.shd_LEstoma_new[j]*LAMBDA(prof.tair_filter[j]+TN0);
		  prof.dLEdz[j][1] = prof.dLAIdz[j] *
		    (solar.prob_beam[j] * (prof.sun_LEstoma[j][1]+prof.sun_LEwet[j][1])
		     + solar.prob_shd[j] * (prof.shd_LEstoma[j][1]+prof.shd_LEwet[j][1]));
		}
	      /* printf("CV48.01 %20.14f %20.14f\n", prof.sun_LEstoma[1][1], prof.sun_LEstoma[jtot][1]); */
              /* printf("CV48.02 %20.14f %20.14f\n", prof.dLEdz[1][1], prof.dLEdz[jtot][1]); */
	    }
	      
	  // compute canopy transpiration and evaporation
	  flux.c_evaporation[1] = 0.;
	  flux.c_transpiration[1] = 0.;
	  temp3 = 0.;
	  for (j=1; j<=jtot; j++)
	    {
	      flux.c_evaporation[1] += prof.dLAIdz[j]
		* (prof.sun_LEwet[j][1]*solar.prob_beam[j]
		   + prof.shd_LEwet[j][1]*solar.prob_shd[j])
		/ LAMBDA(prof.tair_filter[j]+TN0); //convert from W m-2 to kg H2O m-2 s-1
	      flux.c_transpiration[1] += prof.dLAIdz[j] *
		(prof.sun_LEstoma[j][1]*solar.prob_beam[j]
		 + prof.shd_LEstoma[j][1]*solar.prob_shd[j])
		/ LAMBDA(prof.tair_filter[j]+TN0);
	      temp3 += prof.dLAIdz[j]
		* (solar.rnet_sun[j]*solar.prob_beam[j]
		   + solar.rnet_shd[j]*solar.prob_shd[j]);
	    }

	  // total flux
	  flux.c_evapotranspiration[1] = 
	    flux.c_evaporation[1] + flux.c_transpiration[1];
	  flux.evapotranspiration[1] = 
	    flux.c_evapotranspiration[1] + flux.s_evap[1];
	  /* printf("CV48.03 %20.14f %20.14f\n", flux.c_evaporation[1], flux.c_evapotranspiration[1]); */
	      
	  /*
	    Compute air temperature profiles from dispersion matrix

	    Adjust dHdz[1] for soil heat flux

	    sign convention used: fluxes from surface are positive
	    those toward the surface are negative

	    dHdz is net sensible heat exchange per unit leaf area

	    Using filtered temperatures to minimize the system from being mathematically unstable.
	  */

	  // filter temperatures with each interation to minimize numerical instability

	  // Matthias, constant filter
	  //fact.a_filt = 0.85;

	  // Matthias, force conversion with narrowing filter
	  //fact.a_filt = 0.99 - 0.98*((double) i_count)/((double) (i_max-1));

	  // Matthias, narrowing filter to 0.5 at i_max/2
	  fact.a_filt = 0.85 - 0.7*((double) i_count)/((double) (i_max-1));
	  /* printf("CV48.04 %20.14f\n", fact.a_filt); */

	  // Matthias, a_filt=0 from 10 steps before max iteration
	  //   the noise one sees should be filtered off as well
	  //fact.a_filt = __max(0.99 - 1.0*((double) i_count)/((double) (i_max-10)), 0.);

	  // Matthias, original filter
	  // 	      if (i_count < 10)
	  // 	      fact.a_filt = 0.85;
	  //               else
	  // 	      fact.a_filt = 0.5;

	  // conc, for temperature profiles using source/sinks

	  // inputs are source/sink[], scalar_profile[],
	  // ref_val, boundary_flux, unit conversions
	  CONC(prof.dHdz, prof.tair, input.ta, soil.heat, fact.heatcoef);
          /* printf("CV49.01 %20.14f %20.14f\n", prof.dHdz[1], prof.dHdz[jtot]); */
          /* printf("CV49.02 %20.14f %20.14f %20.14f\n", input.ta, soil.heat, fact.heatcoef); */
          /* printf("CV49.03 %20.14f %20.14f\n", prof.tair[1], prof.tair[jtot3-1]); */

	  // filter temperatures to remove numerical instabilities
	  // for each iteration

	  for (j=1; j<=jtot3; j++)
	    if (prof.tair[j] < -30. || prof.tair[j] > 60.)
	      prof.tair[j] = input.ta;
          /* printf("CV49.04 %20.14f %20.14f\n", prof.tair[1], prof.tair[jtot3-1]); */

	  /*
	    Compute vapor density profiles from Dispersion matrix

	    sign convention used: fluxes from surface are positive
	    those toward the surface are negative

	    dLEdZ is net latent heat exchange per unit leaf area
	  */

	  // turbulent transport of H2O

	  for (j=1; j<=jtot; j++)
	    tmp2[j] = prof.dLEdz[j][1]/LAMBDA(prof.tair_filter_save[j]+TN0);
	  CONC(tmp2, tmp1, met.rhova_kg, flux.s_evap[1], 1.);
	  for (j=1; j<=jtot3; j++) prof.rhov_air[j][1] = tmp1[j];
          /* printf("CV49.05 %20.14f %20.14f\n", tmp2[1], tmp2[jtot]); */
          /* printf("CV49.06 %20.14f %20.14f %20.14f\n", met.rhova_kg, flux.s_evap[1], 1.); */
          /* printf("CV49.07 %20.14f %20.14f\n", prof.rhov_air[1][1], prof.rhov_air[jtot3-1][1]); */

	  // filter humidity computations

	  for (j=1; j<=jtot3; j++)
	    {
	      if ((prof.rhov_air[j][1] < 0.) || (prof.rhov_air[j][1] > 0.03))
		prof.rhov_air[j][1] = met.rhova_kg;
	      ea = prof.rhov_air[j][1]/2.165*(prof.tair[j]+TN0);
	      // vapor pressure deficit in [kPa]
	      prof.vpd_air[j] = ES(prof.tair[j]+TN0) - ea;
	    }
          /* printf("CV49.08 %20.14f %20.14f\n", prof.rhov_air[1][1], prof.rhov_air[jtot3-1][1]); */
          /* printf("CV49.09 %20.14f %20.14f\n", prof.vpd_air[1], prof.vpd_air[jtot3-1]); */

	  // Implicit water isotopes diagnostics
	  if (set_switch.wiso == 1 && wiso.implicit == 1)
	    {	      
	      // isotope soil water flux
	      SOIL_FLUX_WISO();
	      /* printf("CV52.01 %20.14f %20.14f\n", flux.soilevap[1], flux.soilevap[wiso.nwater]); */
	      /* printf("CV52.02 %20.14f %20.14f\n", flux.litterevap[1], flux.litterevap[wiso.nwater]); */
	      /* printf("CV52.03 %20.14f %20.14f\n", flux.s_evap[1], flux.s_evap[wiso.nwater]); */
	      
	      // leaf water enrichment
	      LEAF_WISO();
	      /* printf("CV52.04 %20.14f %20.14f\n", prof.sun_wi[1], prof.sun_wi[jtot]); */
	      /* printf("CV52.05 %20.14f %20.14f\n", prof.shd_wi[1], prof.shd_wi[jtot]); */
	      /* printf("CV52.06 %20.14f %20.14f\n", prof.wa[1], prof.wa[jtot]); */
	      /* printf("CV52.07 %20.14f %20.14f\n", prof.sun_h[1], prof.sun_h[jtot]); */
	      /* printf("CV52.08 %20.14f %20.14f\n", prof.rs_fact[1], prof.rs_fact[jtot]); */
	      /* printf("CV52.09 %20.14f %20.14f\n", prof.sun_gross[1], prof.sun_gross[jtot]); */
	      /* printf("CV52.10 %20.14f %20.14f\n", prof.shd_gross[1], prof.shd_gross[jtot]); */
	      /* printf("CV52.11 %20.14f %20.14f\n", prof.sun_LEstoma_new[1], prof.sun_LEstoma_new[jtot]); */
	      /* printf("CV52.12 %20.14f %20.14f\n", prof.shd_LEstoma_new[1], prof.shd_LEstoma_new[jtot]); */

	      // isotope canopy transpiration and evaporation
	      CANOPY_FLUX_WISO();
              /* printf("CV52.13 %20.14f %20.14f\n", prof.sun_alpha_k[1][1], prof.sun_alpha_k[jtot][wiso.nwater]); */
              /* printf("CV52.14 %20.14f %20.14f\n", prof.shd_alpha_k[1][1], prof.shd_alpha_k[jtot][wiso.nwater]); */
              /* printf("CV52.15 %20.14f %20.14f\n", prof.sun_alpha_equ[1][1], prof.sun_alpha_equ[jtot][wiso.nwater]); */
              /* printf("CV52.16 %20.14f %20.14f\n", prof.shd_alpha_equ[1][1], prof.shd_alpha_equ[jtot][wiso.nwater]); */
              /* printf("CV52.17 %20.14f %20.14f\n", prof.sun_peclet[1][1], prof.sun_peclet[jtot][wiso.nwater]); */
              /* printf("CV52.18 %20.14f %20.14f\n", prof.shd_peclet[1][1], prof.shd_peclet[jtot][wiso.nwater]); */
              /* printf("CV52.19 %20.14f %20.14f\n", prof.sun_fem[1][1], prof.sun_fem[jtot][wiso.nwater]); */
              /* printf("CV52.20 %20.14f %20.14f\n", prof.shd_fem[1][1], prof.shd_fem[jtot][wiso.nwater]); */
              /* printf("CV52.21 %20.14f %20.14f\n", prof.sun_craig[1][1], prof.sun_craig[jtot][wiso.nwater]); */
              /* printf("CV52.22 %20.14f %20.14f\n", prof.shd_craig[1][1], prof.shd_craig[jtot][wiso.nwater]); */
              /* printf("CV52.23 %20.14f %20.14f\n", prof.sun_leafwater_e[1][1], prof.sun_leafwater_e[jtot][wiso.nwater]); */
              /* printf("CV52.24 %20.14f %20.14f\n", prof.shd_leafwater_e[1][1], prof.shd_leafwater_e[jtot][wiso.nwater]); */
              /* printf("CV52.25 %20.14f %20.14f\n", prof.sun_leafwater[1][1], prof.sun_leafwater[jtot][wiso.nwater]); */
              /* printf("CV52.26 %20.14f %20.14f\n", prof.shd_leafwater[1][1], prof.shd_leafwater[jtot][wiso.nwater]); */
              /* printf("CV52.27 %20.14f %20.14f\n", prof.sun_trans_rtrans[1][1], prof.sun_trans_rtrans[jtot][wiso.nwater]); */
              /* printf("CV52.28 %20.14f %20.14f\n", prof.shd_trans_rtrans[1][1], prof.shd_trans_rtrans[jtot][wiso.nwater]); */
              /* printf("CV52.29 %20.14f %20.14f\n", prof.sun_rtrans[1][1], prof.sun_rtrans[jtot][wiso.nwater]); */
              /* printf("CV52.30 %20.14f %20.14f\n", prof.shd_rtrans[1][1], prof.shd_rtrans[jtot][wiso.nwater]); */
    	      
	      // turbulent transport of H2O18, DHO and H216O
	      for (mc=2; mc<=wiso.nwater+1; mc++)
		{
		  // Vapour input in layers
		  for (j=1; j<=jtot; j++)
		    {
		      tmp2[j] = prof.dLEdz[j][mc]/LAMBDA(prof.tair_filter_save[j]+TN0);
#ifdef DIAG
		      if (mc == 4 && abs(DELTA1000_H2O(prof.dLEdz[j][mc], prof.dLEdz[j][1], mc)) >= isotolerance)
			printf("Trans happend at layer %i during time step %ld: %g\n",
			       j, time_var.daytime,
			       DELTA1000_H2O(prof.dLEdz[j][mc], prof.dLEdz[j][1], mc));
#endif
		    }
		  // Vapour input at top
		  tmp3 = ALPHA_EQU_H2O(input.ta+TN0, mc)
		    * INVDELTA1000_H2O(input.dppt[mc], mc) * met.rhova_kg; // in equi with rain
#ifdef DIAG
		  if (mc == 4 && abs(DELTA1000_H2O(tmp3, met.rhova_kg, mc)) >= isotolerance)
		    printf("Input happend during time step %ld: %g\n",
			   time_var.daytime, DELTA1000_H2O(tmp3, met.rhova_kg, mc));
		  if (mc == 4 && abs(DELTA1000_H2O(flux.s_evap[mc], flux.s_evap[1], mc)) >= isotolerance)
		    printf("Evap happend during time step %ld: %g\n",
			   time_var.daytime, DELTA1000_H2O(flux.s_evap[mc], flux.s_evap[1], mc));
#endif
		  CONC(tmp2, tmp1, tmp3, flux.s_evap[mc], 1.);
		  for (j=1; j<=jtot3; j++) 
		    prof.rhov_air[j][mc] = tmp1[j];
		}
	      
	      // filter humidity computations
	      for (j=1; j<=jtot3; j++)
		if (prof.rhov_air[j][1] == met.rhova_kg && j != jtot3)
		  for (mc=2; mc<=wiso.nwater+1; mc++)
		    prof.rhov_air[j][mc] = prof.rhov_air[j][1] * prof.rvapour[j][mc];
	    } // implicit water isotope diagnostics
	      
              /*
                sign convention used: photosynthetic uptake is positive
                respiration is negative

                prof.dPsdz is net photosynthesis per unit leaf area,
                the profile was converted to units of mg m-2 s-1 to be
                consistent with inputs to CONC
              */

              // change sign of dPsdz
	  for (j=1; j<=jtot; j++)
	    prof.source_co2[j] = -prof.dPsdz[j];
          /* printf("CV50.01 %20.14f %20.14f\n", prof.source_co2[1], prof.source_co2[jtot]); */

	  // compute bole respiration
	  // It is added to the prof.source_CO2 array.
	  BOLE_RESPIRATION();
          /* printf("CV50.02 %20.14f %20.14f\n", bole.calc, bole.factor); */
          /* printf("CV50.03 %20.14f %20.14f\n", bole.respiration_mole, bole.respiration_mg); */
          /* printf("CV50.04 %20.14f %20.14f\n", bole.layer[1], bole.layer[jtot]); */
          /* printf("CV50.05 %20.14f %20.14f\n", prof.source_co2[1], prof.source_co2[jtot]); */
              
	  // compute soil respiration
	  SOIL_RESPIRATION();
          /* printf("CV50.06 %20.14f %20.14f\n", soil.respiration_mole, soil.respiration_mg); */
          /* printf("CV50.07 %20.14f %20.14f\n", soil.respiration_auto, soil.respiration_hetero); */

	  // to convert umol m-3 to umol/mol we have to consider
	  // Pc/Pa = [CO2]ppm = rhoc ma/ rhoa mc
	  fact.co2 = (mass_air/mass_CO2)*met.air_density_mole;

	  CONC(prof.source_co2, prof.co2_air, input.co2air, soil.respiration_mole, fact.co2);
          /* printf("CV50.08 %20.14f %20.14f\n", prof.source_co2[1], prof.source_co2[jtot]); */
          /* printf("CV50.09 %20.14f %20.14f %20.14f\n", input.co2air, soil.respiration_mole, fact.co2); */
          /* printf("CV50.10 %20.14f %20.14f\n", prof.co2_air[1], prof.co2_air[jtot3-1]); */

	  // Integrate source-sink strengths to estimate canopy flux

	  sumh = 0.; // sensible heat
	  sumle = 0.; // latent heat
	  sumrn = 0.; // net radiation
	  sumlout = 0.; // radiation
	  can_ps_mol = 0.; // canopy photosynthesis
	  can_gpp = 0.; // canopy GPP = can_ps_mol + dark respiration
	  canresp = 0.; // canopy respiration
	  sumlai = 0.; // leaf area
	  sumksi = 0.; // canopy stomatal conductance
	  tleaf_mean = 0.; // mean leaf temperature
	  tavg_sun = 0.; // avg sunlit temperature
	  tavg_shd = 0.; // avg shaded temperature
	  sumbolelay = 0.; // integrated bole respiration
	  sum_cica = 0.;
	  sum_gs = 0.;

	  for (j=1; j<=jtot; j++)
	    {
	      sumh += prof.dHdz[j];
	      sumle += prof.dLEdz[j][1];
	      sumrn += prof.dRNdz[j];
	      sumlout += prof.dLoutdz[j];
	      can_ps_mol += prof.dPsdz[j];
	      can_gpp += prof.dGPPdz[j];
	      canresp += prof.dRESPdz[j];
	      sumksi += prof.dStomCondz[j];
	      sumlai += prof.dLAIdz[j];
	      tleaf_mean += prof.sun_tleaf[j]*solar.prob_beam[j] + prof.shd_tleaf[j]*solar.prob_shd[j];
	      prof.tleaf[j] = prof.sun_tleaf[j]*solar.prob_beam[j] + prof.shd_tleaf[j]*solar.prob_shd[j]; //Tleaf per layer (sun and shade)
	      sumbolelay += bole.layer[j];

	      // need to weight by sun and shaded leaf areas then divide by LAI

	      tavg_sun += prof.sun_tleaf[j]*prof.dLAIdz[j];
	      tavg_shd += prof.shd_tleaf[j]*prof.dLAIdz[j];
	    }
          /* printf("CV50.11 %20.14f %20.14f %20.14f\n", sumh, sumle, sumrn); */
          /* printf("CV50.12 %20.14f %20.14f %20.14f\n", sumlout, can_ps_mol, can_gpp); */
          /* printf("CV50.13 %20.14f %20.14f %20.14f\n", canresp, sumksi, sumlai); */
	  /* printf("T: %i %20.14f %20.14f\n", i_count, prof.shd_tleaf[1], prof.shd_tleaf[40]); */
          /* printf("CV50.14 %20.14f %20.14f\n", prof.tleaf[1], prof.tleaf[jtot]); */
          /* printf("CV50.15 %20.14f %20.14f\n", tavg_sun, tavg_shd); */
	  ebalance = sumrn - sumle - sumh;

	  flux.photosyn = can_ps_mol;
          /* printf("CV50.16 %20.14f %20.14f\n", ebalance, flux.photosyn); */

	  // calculate gpp : can_ps_mol = photosynthesis - photorespiration - dark respiration or can_ps_mol = GPP - leaf respiration
	  //can_gpp = can_ps_mol + canresp;

	  // mean canopy leaf temperature

	  tleaf_mean /= jtot;

	  // leaf area weighted temperatures

	  tavg_sun /= sumlai;
	  tavg_shd /= sumlai;
	  /* printf("CV50.017 %20.14f %20.14f %20.14f\n", tleaf_mean, tavg_sun, tavg_shd); */

	  // Energy exchanges at the soil

	  sbalance = soil.rnet - soil.lout
	    - soil.soilevap - soil.litterevap
	    - soil.heat - soil.gsoil;

	  rnet_soil = soil.rnet - soil.lout;
          /* printf("CV50.18 %20.14f %20.14f\n", rnet_soil, sbalance); */

	  // canopy scale flux densities, vegetation plus soil


	  sumh += soil.heat;
	  sumle += soil.evap;
	  sumrn += rnet_soil;
	  sumlout += soil.lout;
	  temp3 += soil.rnet;
	      
	  tbalance = sbalance + ebalance;

	  met.H = sumh;
          /* printf("CV50.19 %20.14f %20.14f %20.14f\n", sumh, sumle, sumrn); */
          /* printf("CV50.20 %20.14f %20.14f %20.14f\n", sumlout, temp3, tbalance); */
          /* printf("CV50.21 %20.14f\n", met.H); */

	  // filter iterative variables

	  // soil temp and fluxes
	  soil.tsfc_filter       = fact.a_filt * soil.tsfc
	    + (1.-fact.a_filt) * soil.tsfc_filter;
	  soil.soilevap_filter   = fact.a_filt * soil.soilevap
	    + (1.-fact.a_filt) * soil.soilevap_filter;
	  soil.litterevap_filter = fact.a_filt * soil.litterevap
	    + (1.-fact.a_filt) * soil.litterevap_filter;
	  soil.T_l_filter       = fact.a_filt * soil.T_l
	    + (1.-fact.a_filt) * soil.T_l_filter;
	  for (j=0; j<=soil.n_soil+1; j++)
	    soil.T_soil_filter[j] = fact.a_filt * soil.T_soil[j]
	      + (1.-fact.a_filt) * soil.T_soil_filter[j];
	  /* printf("CV50.22 %20.14f %20.14f %20.14f\n", soil.tsfc_filter, soil.soilevap_filter, soil.litterevap_filter); */
          /* printf("CV50.23 %20.14f\n", soil.T_l_filter); */
          /* printf("CV50.231 %20.14f %20.14f\n", soil.T_soil_filter[0], soil.T_soil_filter[soil.n_soil+1]); */
	  // met variables
	  met.H_filter           = fact.a_filt * met.H
	    + (1.-fact.a_filt) * met.H_filter;
	  met.ustar_filter       = fact.a_filt * met.ustar
	    + (1.-fact.a_filt) * met.ustar_filter;
          /* printf("CV50.24 %20.14f %20.14f\n", met.H_filter, met.ustar_filter); */
	  // air variables
	  for (j=1; j<=jtot3; j++)
	    {
	      prof.tair_filter[j] = fact.a_filt * prof.tair[j]
		+ (1.-fact.a_filt) * prof.tair_filter[j];
	      prof.rhov_air_filter[j][1] = fact.a_filt * prof.rhov_air[j][1]
		+ (1.-fact.a_filt) * prof.rhov_air_filter[j][1];
	      prof.co2_air_filter[j] = fact.a_filt * prof.co2_air[j]
		+ (1.-fact.a_filt) * prof.co2_air_filter[j];
	    }
          /* printf("CV50.25 %20.14f %20.14f\n", prof.tair_filter[1], prof.tair_filter[jtot3-1]); */
          /* printf("CV50.26 %20.14f %20.14f\n", prof.rhov_air_filter[1][1], prof.rhov_air_filter[jtot3-1][1]); */
          /* printf("CV50.27 %20.14f %20.14f\n", prof.co2_air_filter[1], prof.co2_air_filter[jtot3-1]); */
	  // leaf variables
	  for (j=1; j<=jtot; j++)
	    {
	      prof.sun_tleaf_filter[j] = fact.a_filt * prof.sun_tleaf[j]
		+ (1.-fact.a_filt) * prof.sun_tleaf_filter[j];
	      prof.shd_tleaf_filter[j] = fact.a_filt * prof.shd_tleaf[j]
		+ (1.-fact.a_filt) * prof.shd_tleaf_filter[j];
	      prof.wet_coef_filter[j] = fact.a_filt * prof.wet_coef[j]
		+ (1.-fact.a_filt) * prof.wet_coef_filter[j];
	    }
          /* printf("CV50.28 %20.14f %20.14f\n", prof.sun_tleaf_filter[1], prof.sun_tleaf_filter[jtot]); */
          /* printf("CV50.29 %20.14f %20.14f\n", prof.shd_tleaf_filter[1], prof.shd_tleaf_filter[jtot]); */
          /* printf("CV50.30 %20.14f %20.14f\n", prof.wet_coef_filter[1], prof.wet_coef_filter[jtot]); */

	  // Implicit water isotopes diagnostics
	  if (set_switch.wiso == 1 && wiso.implicit == 1)
	    {	      
	      for (mc=2; mc<=wiso.nwater+1; mc++)
		for (j=1; j<=jtot3; j++)
		  prof.rhov_air_filter[j][mc] = fact.a_filt*prof.rhov_air[j][mc]
		    + (1.-fact.a_filt)*prof.rhov_air_filter[j][mc];
		  
	      // isotopes in vapour
	      for (mc=2; mc<=wiso.nwater+1; mc++)
		for (j=1; j<=jtot3; j++)
		  {
#ifdef DIAG
		    if (fabs(prof.rhov_air_filter[j][1]) < 1e-15 &&
			(prof.rhov_air_filter[j][mc] != 0. || prof.rhov_air_filter[j][1] != 0.))
		      printf("\nBefore Isorat 1 tracer %i @ %ld: rare %g abundant %g\n", 
			     mc, time_var.daytime, 
			     prof.rhov_air_filter[j][mc], prof.rhov_air_filter[j][1]);
#endif
		    prof.rvapour[j][mc] = ISORAT(&prof.rhov_air_filter[j][mc], &prof.rhov_air_filter[j][1],
							&wiso.lost[mc], &wiso.lost[1], mc);
		    if (wiso.lost[mc] != 0. && soil.lost0 == 0)
		      {
			printf("\nLost 04.1 @ tracer %i t %ld\n", mc, time_var.daytime);
			soil.lost0 = 1;
		      }
		  }
	    }

	  // Net radiation balance at the top of the canopy

	  netrad = (solar.beam_flux_par[jtot1] + solar.par_down[jtot1] - solar.par_up[jtot1]) / 4.6 // PAR
	    + solar.beam_flux_nir[jtot1] + solar.nir_dn[jtot1] - solar.nir_up[jtot1] // NIR
	    + solar.ir_dn[jtot1] - solar.ir_up[jtot1]; // IR
          /* printf("CV50.31 %20.14f\n", netrad); */

	  // test for convergenece between the sum of the net radiation flux profile and the
	  // net flux exiting the canopy

	  //etest = fabs((sumrn-netrad)/sumrn);
	  etest       = (sumrn-netrad)/sumrn;
	  etest_diff1 = etest_old1 - etest;
	  etest_diff2 = etest_old2 - etest_old1;
	  etest_old2  = etest_old1;
	  etest_old1  = etest;

	  itest = 0.;
	  if (set_switch.wiso == 1 && wiso.implicit == 1)
	    itest = DELTA1000_H2O(flux.evapotranspiration[2], 
				  flux.evapotranspiration[1], 2);
	  itest_diff1 = itest_old1 - itest;
	  itest_diff2 = itest_old2 - itest_old1;
	  itest_old2  = itest_old1;
	  itest_old1  = itest;
          /* printf("CV50.32 %20.14f %20.14f %20.14f\n", etest, etest_diff1, etest_diff2); */
          /* printf("CV50.33 %20.14f %20.14f\n", etest_old2, etest_old1); */
          /* printf("CV50.34 %20.14f %20.14f %20.14f\n", itest, itest_diff1, itest_diff2); */
          /* printf("CV50.35 %20.14f %20.14f\n", itest_old2, itest_old1); */

	  i_count++;
	  time_var.count = i_count;
          /* printf("CV50.36 %d\n", time_var.count); */
	}
      //while(etest > .005 && i_count < i_max); // end of while for integrated looping
      while((fabs(etest_diff1) > 0.005 || fabs(etest_diff2) > 0.005
	     || fabs(itest_diff1) > 0.03 || fabs(itest_diff2) > 0.03)
	    && i_count < i_max); // end of while for integrated looping
      //           while((fabs(etest_diff1) > 1e-15 || fabs(etest_diff2) > 1e-15
      // 		 || fabs(itest_diff1) > 1e-15 || fabs(itest_diff2) > 1e-15)
      // 		&& i_count < i_max); // end of while for integrated looping
	  
      totalcount++;
      if (solar.sine_beta <= 0.05)
	ntotalcount++;
      else
	dtotalcount++;
      if (i_count == i_max)
	{
	  count40++;
	  if (solar.sine_beta <= 0.05)
	    ncount40++;
	  else
	    dcount40++;
	}
      if (print_balance == 1)
      	printf("Day %3i Time %5.1f i_c %3i: diff %6.2f%% net %7.2f sum %7.2f, soil %9.2e leaf %9.2e total %9.2e\n",
      	       time_var.days, time_var.local_time, i_count,
      	       etest*100, netrad, sumrn,
      	       sbalance, ebalance, tbalance);

      // Explicit water isotope diagnostics
      if (set_switch.wiso == 1 && wiso.implicit == 0)
	{
	  // isotopes in vapour
	  // like that, we do not have to transport isotopes all the time
	  //   but take the isotopes from one time step before
	  for (mc=2; mc<=wiso.nwater+1; mc++)
	    for (j=1; j<=jtot3; j++)
	      prof.rhov_air_filter[j][mc] = prof.rhov_air_filter_save[j] * prof.rvapour[j][mc];
	      
	  // isotope soil water flux
	  SOIL_FLUX_WISO();
          /* printf("CV51.01 %20.14f %20.14f\n", flux.soilevap[1], flux.soilevap[wiso.nwater]); */
          /* printf("CV51.02 %20.14f %20.14f\n", flux.litterevap[1], flux.litterevap[wiso.nwater]); */
          /* printf("CV51.03 %20.14f %20.14f\n", flux.s_evap[1], flux.s_evap[wiso.nwater]); */
      
	  // leaf water enrichment
	  LEAF_WISO();
          /* printf("CV51.04 %20.14f %20.14f\n", prof.sun_wi[1], prof.sun_wi[jtot]); */
          /* printf("CV51.05 %20.14f %20.14f\n", prof.shd_wi[1], prof.shd_wi[jtot]); */
          /* printf("CV51.06 %20.14f %20.14f\n", prof.wa[1], prof.wa[jtot]); */
          /* printf("CV51.07 %20.14f %20.14f\n", prof.sun_h[1], prof.sun_h[jtot]); */
          /* printf("CV51.08 %20.14f %20.14f\n", prof.rs_fact[1], prof.rs_fact[jtot]); */
          /* printf("CV51.09 %20.14f %20.14f\n", prof.sun_gross[1], prof.sun_gross[jtot]); */
          /* printf("CV51.10 %20.14f %20.14f\n", prof.shd_gross[1], prof.shd_gross[jtot]); */
          /* printf("CV51.11 %20.14f %20.14f\n", prof.sun_LEstoma_new[1], prof.sun_LEstoma_new[jtot]); */
          /* printf("CV51.12 %20.14f %20.14f\n", prof.shd_LEstoma_new[1], prof.shd_LEstoma_new[jtot]); */
  	      
	  // isotope canopy transpiration and evaporation
	  CANOPY_FLUX_WISO();
          /* printf("CV51.13 %20.14f %20.14f\n", prof.sun_alpha_k[1][1], prof.sun_alpha_k[jtot][wiso.nwater]); */
          /* printf("CV51.14 %20.14f %20.14f\n", prof.shd_alpha_k[1][1], prof.shd_alpha_k[jtot][wiso.nwater]); */
          /* printf("CV51.15 %20.14f %20.14f\n", prof.sun_alpha_equ[1][1], prof.sun_alpha_equ[jtot][wiso.nwater]); */
          /* printf("CV51.16 %20.14f %20.14f\n", prof.shd_alpha_equ[1][1], prof.shd_alpha_equ[jtot][wiso.nwater]); */
          /* printf("CV51.17 %20.14f %20.14f\n", prof.sun_peclet[1][1], prof.sun_peclet[jtot][wiso.nwater]); */
          /* printf("CV51.18 %20.14f %20.14f\n", prof.shd_peclet[1][1], prof.shd_peclet[jtot][wiso.nwater]); */
          /* printf("CV51.19 %20.14f %20.14f\n", prof.sun_fem[1][1], prof.sun_fem[jtot][wiso.nwater]); */
          /* printf("CV51.20 %20.14f %20.14f\n", prof.shd_fem[1][1], prof.shd_fem[jtot][wiso.nwater]); */
          /* printf("CV51.21 %20.14f %20.14f\n", prof.sun_craig[1][1], prof.sun_craig[jtot][wiso.nwater]); */
          /* printf("CV51.22 %20.14f %20.14f\n", prof.shd_craig[1][1], prof.shd_craig[jtot][wiso.nwater]); */
          /* printf("CV51.23 %20.14f %20.14f\n", prof.sun_leafwater_e[1][1], prof.sun_leafwater_e[jtot][wiso.nwater]); */
          /* printf("CV51.24 %20.14f %20.14f\n", prof.shd_leafwater_e[1][1], prof.shd_leafwater_e[jtot][wiso.nwater]); */
          /* printf("CV51.25 %20.14f %20.14f\n", prof.sun_leafwater[1][1], prof.sun_leafwater[jtot][wiso.nwater]); */
          /* printf("CV51.26 %20.14f %20.14f\n", prof.shd_leafwater[1][1], prof.shd_leafwater[jtot][wiso.nwater]); */
          /* printf("CV51.27 %20.14f %20.14f\n", prof.sun_trans_rtrans[1][1], prof.sun_trans_rtrans[jtot][wiso.nwater]); */
          /* printf("CV51.28 %20.14f %20.14f\n", prof.shd_trans_rtrans[1][1], prof.shd_trans_rtrans[jtot][wiso.nwater]); */
          /* printf("CV51.29 %20.14f %20.14f\n", prof.sun_rtrans[1][1], prof.sun_rtrans[jtot][wiso.nwater]); */
          /* printf("CV51.30 %20.14f %20.14f\n", prof.shd_rtrans[1][1], prof.shd_rtrans[jtot][wiso.nwater]); */
  	      
	  // turbulent transport of H2O18, DHO and H216O
	  for (mc=2; mc<=wiso.nwater+1; mc++)
	    {
	      // Vapour input in layers
	      for (j=1; j<=jtot; j++)
		{
		  tmp2[j] = prof.dLEdz[j][mc]/LAMBDA(prof.tair_filter_save[j]+TN0);
#ifdef DIAG
		  if (mc == 4 && abs(DELTA1000_H2O(prof.dLEdz[j][mc], prof.dLEdz[j][1], mc)) >= isotolerance)
		    printf("Trans happend at layer %i during time step %ld: %g\n",
			   j, time_var.daytime,
			   DELTA1000_H2O(prof.dLEdz[j][mc], prof.dLEdz[j][1], mc));
#endif
		}
	      // Vapour input at top
	      tmp3 = ALPHA_EQU_H2O(prof.tair_filter_save[j]+TN0, mc)
		* INVDELTA1000_H2O(input.dppt[mc], mc) * met.rhova_kg; // in equi with rain
#ifdef DIAG
	      if (mc == 4 && abs(DELTA1000_H2O(tmp3, met.rhova_kg, mc)) >= isotolerance)
		printf("Input happend during time step %ld: %g\n",
		       time_var.daytime, DELTA1000_H2O(tmp3, met.rhova_kg, mc));
	      if (mc == 4 && abs(DELTA1000_H2O(flux.s_evap[mc], flux.s_evap[1], mc)) >= isotolerance)
		printf("Evap happend during time step %ld: %g\n",
		       time_var.daytime, DELTA1000_H2O(flux.s_evap[mc], flux.s_evap[1], mc));
#endif
	      CONC(tmp2, tmp1, tmp3, flux.s_evap[mc], 1.);
	      for (j=1; j<=jtot3; j++) 
		prof.rhov_air[j][mc] = tmp1[j];
	    }
	      
	  // filter humidity computations
	  for (j=1; j<=jtot3; j++)
	    if (prof.rhov_air[j][1] == met.rhova_kg && j != jtot3)
	      for (mc=2; mc<=wiso.nwater+1; mc++)
		prof.rhov_air[j][mc] = prof.rhov_air[j][1] * prof.rvapour[j][mc];

	  // isotopes in vapour
	  for (mc=2; mc<=wiso.nwater+1; mc++)
	    for (j=1; j<=jtot3; j++)
	      { // take rhov_air not rhov_air_filter
#ifdef DIAG
		if (fabs(prof.rhov_air[j][1]) < 1e-15 &&
		    (prof.rhov_air[j][mc] != 0. || prof.rhov_air[j][1] != 0.))
		  printf("\nBefore Isorat 2 tracer %i @ %ld: rare %g abundant %g\n", 
			 mc, time_var.daytime,
			 prof.rhov_air[j][mc], prof.rhov_air[j][1]);
#endif
		prof.rvapour[j][mc] = ISORAT(&prof.rhov_air[j][mc], &prof.rhov_air[j][1],
						    &wiso.lost[mc], &wiso.lost[1], mc);
		if (wiso.lost[mc] != 0. && soil.lost0 == 0)
		  {
		    printf("\nLost 04.1 @ tracer %i t %ld\n", mc, time_var.daytime);
		    soil.lost0 = 1;
		  }
	      }
	} // explicit water isotope diagnostics

      // subtract evaporation from wet leaf surfaces from canopy water storage
      temp1=0;
      for (j=1; j<=jtot; j++)
	{
	  if (prof.cws[j][1] > 0.)
	    {
	      for (mc=2; mc<=wiso.nwater+1; mc++)
		{
#ifdef DIAG
		  if (prof.cws[j][1] != 0. || prof.cws[j][1] != 0.)
		    if (fabs(prof.cws[j][mc]) < 1e-15 &&
			(prof.cws[j][mc] != 0. || prof.cws[j][1] != 0.))
		      printf("\nBefore Isorat 3 tracer %i @ %ld: rare %g abundant %g\n", 
			     mc, time_var.daytime, 
			     prof.cws[j][mc], prof.cws[j][1]);
#endif
		  rcws[mc] = ISORAT(&prof.cws[j][mc], &prof.cws[j][1],
				    &wiso.lost[mc], &wiso.lost[1], mc);
		  if (wiso.lost[mc] != 0. && soil.lost0 == 0)
		    {
		      printf("\nLost 01 @ tracer %i t %ld\n",mc,time_var.daytime);
		      soil.lost0 = 1;
		    }
		}
	      prof.cws[j][1] = __max(0., prof.cws[j][1] 
					    - (prof.sun_LEwet[j][1]*solar.prob_beam[j]
					       + prof.shd_LEwet[j][1]*solar.prob_shd[j])
					    * prof.dLAIdz[j]
					    / LAMBDA(prof.tair_filter_save[j] + TN0) * time_var.time_step);
	      // No fractionation for interception evaporation
	      for (mc=2; mc<=wiso.nwater+1; mc++)
		prof.cws[j][mc] = rcws[mc] * prof.cws[j][1];
	      // Check: if this happens, we have to code some thresholds here
	      for (mc=1; mc<=wiso.nwater+1; mc++)
		if (prof.cws[j][mc] < 0.)
		  printf("\nC1<0 trac %i layer %i: %e\n", mc, j, prof.cws[j][mc]);
	    }
	  temp1 += prof.cws[j][1];
	}

      c_transpiration_mole = 1000.*1000.* flux.c_transpiration[1]/18.; // mmol m-2 s-1

      if (soil.camillo == 0)
	{
	  // compute litter moisture and litter drainage
	  LITTER_H2O();
          /* printf("CV53.01 %20.14f %20.14f\n", flux.soilinfl[1], flux.soilinfl[wiso.nwater]); */
          /* printf("CV53.02 %20.14f %20.14f\n", soil.qinfl[1], soil.qinfl[wiso.nwater]); */
          /* printf("CV53.03 %20.14f %20.14f\n", soil.theta_l[1], soil.theta_l[wiso.nwater]); */
	  // compute soil moisture in different layers
	  SOIL_H2O();
          /* printf("CV54.01 %20.14f %20.14f\n", soil.qseva[1], soil.qseva[wiso.nwater]); */
          /* printf("CV54.03 %20.14f %20.14f\n", soil.qdrai[1], soil.qdrai[wiso.nwater]); */
          /* printf("CV54.04 %20.14f %20.14f\n", soil.swp[1], soil.swp[soil.n_soil]); */
          /* printf("CV54.05 %20.14f %20.14f\n", soil.swp_mm[1], soil.swp_mm[soil.n_soil]); */
          /* printf("CV54.06 %20.14f %20.14f\n", soil.k_theta[1], soil.k_theta[soil.n_soil]); */
          /* printf("CV54.07 %20.14f %20.14f\n", soil.a[1], soil.a[soil.n_soil]); */
          /* printf("CV54.08 %20.14f %20.14f\n", soil.b[1], soil.b[soil.n_soil]); */
          /* printf("CV54.09 %20.14f %20.14f\n", soil.c[1], soil.c[soil.n_soil]); */
          /* printf("CV54.10 %20.14f %20.14f\n", soil.r[1][1], soil.r[soil.n_soil][wiso.nwater]); */
          /* printf("CV54.11 %20.14f %20.14f %20.14f\n", soil.soil_mm, soil.soil_mm_root, soil.soil_mm_50); */
          /* printf("CV54.13 %20.14f %20.14f\n", soil.theta[1][1], soil.theta[soil.n_soil][wiso.nwater]); */

	}

      // check and convert units of all components to desirable values
      // Convert to umol m-2 s-1
      can_ps_mg = can_ps_mol * mass_CO2/ 1000.;

      wue = can_ps_mg / (1000. * flux.c_transpiration[1]); /* mg co2 g h20 */

      fc_mg = -(can_ps_mg - soil.respiration_mg - bole.respiration_mg);
      fc_mol = fc_mg * 1000./mass_CO2;

      sumA     = 0.;
      sum_cica = 0.;
      sum_gs   = 0.;
      for (j=1; j<=jtot; j++)
	{
	  sun_A     = __max(prof.dLAIdz[j]*prof.sun_A[j]*solar.prob_beam[j], 0.);
	  shd_A     = __max(prof.dLAIdz[j]*prof.shd_A[j]*solar.prob_shd[j], 0.);
	  sumA     += sun_A + shd_A;
	  sum_cica += sun_A*prof.sun_cica[j] + shd_A*prof.shd_cica[j];
	  sum_gs   += sun_A*prof.sun_gs_mol[j] + shd_A*prof.shd_gs_mol[j];
	}
      if (sumA > 0.)
	{
	  ave_cica = sum_cica/sumA;
	  ave_gs = sum_gs/sumA;
	}
      else
	{
	  ave_cica = undef;
	  ave_gs = undef;
	}

      if (set_switch.d13c == 1)
	{
	  // call carbon isotope subroutine
	  // it produces outputs based on inputs of Ci/Ca, gs, A, Tl
	  CARBON_ISOTOPE();

	  sun_A  = 0.;
	  shd_A  = 0.;
	  sumA   = 0.;
	  sum13C = 0.; // integrated source fluxe of 13C
	  sumD13 = 0.; // avg 13C discrimination, D13, weighted by photosynthesis
	  sumD13_long = 0.;
	  sumD13_long_isotope = 0.;

	  ave_D13_a = 0;
	  sumD13_a = 0;
	  ave_D13_ab = 0;
	  sumD13_ab = 0;
	  ave_D13_asal = 0;
	  sumD13_asal = 0;
	  ave_D13_b = 0;
	  sumD13_b = 0;
	  ave_csca = 0;
	  sum_csca=0;
	  ave_cica = 0;
	  sum_cica=0;
	  ave_ccca = 0;
	  sum_ccca=0;
	  ave_gs = 0;
	  sum_gs=0;
	  ave_gm = 0;
	  sum_gm=0;

	  for (j=1; j<=jtot; j++)
	    sum13C += prof.sour13co2[j];

	  for (j=1; j<=jtot; j++)
	    {
	      sun_A = __max(prof.dLAIdz[j]*prof.sun_A[j]*solar.prob_beam[j], 0.);
	      shd_A = __max(prof.dLAIdz[j]*prof.shd_A[j]*solar.prob_shd[j], 0.);
	      sumA += sun_A + shd_A;
	      sumD13 += sun_A*prof.sun_D13[j] + shd_A*prof.shd_D13[j];
	      sumD13_long += sun_A*prof.sun_D13_long[j] + shd_A*prof.shd_D13_long[j];
	      sumD13 += sun_A*prof.sun_D13[j] +	shd_A*prof.shd_D13[j];
	      sumD13_long += sun_A*prof.sun_D13_long[j] + shd_A*prof.shd_D13_long[j];
	      sumD13_a += sun_A*prof.sun_D13_a[j] + shd_A*prof.shd_D13_a[j];
	      sumD13_ab += sun_A*prof.sun_D13_ab[j] + shd_A*prof.shd_D13_ab[j];
	      sumD13_asal += sun_A*prof.sun_D13_asal[j] + shd_A*prof.shd_D13_asal[j];
	      sumD13_b += sun_A*prof.sun_D13_b[j] + shd_A*prof.shd_D13_b[j];
	      sumD13_long_isotope += prof.D13C_long[j]*prof.dPsdz[j]; // ???

	      sum_csca += (sun_A + shd_A)*prof.csca[j];
	      sum_cica += sun_A*prof.sun_cica[j] + shd_A*prof.shd_cica[j];
	      sum_ccca += (sun_A + shd_A)*prof.ccca[j];
	      sum_gs += sun_A*prof.sun_gs_mol[j] + shd_A*prof.shd_gs_mol[j];
	      //MC20110804 sum_gm += sumA*prof.gm[j];
	      sum_gm += (sun_A+shd_A)*prof.gm[j];
	    }

	  // Del 13 weighted by Photosynthesis

	  if (sumA > 0.)
	    {
	      ave_D13 = sumD13/sumA;
	      ave_D13_long = sumD13_long/sumA;
	      ave_D13_a = sumD13_a/sumA;
	      ave_D13_ab = sumD13_ab/sumA;
	      ave_D13_asal = sumD13_asal/sumA;
	      ave_D13_b = sumD13_b/sumA;
		  
	      ave_D13_long_isotope = sumD13_long_isotope/sumA;

	      ave_csca = sum_csca/sumA;
	      ave_cica = sum_cica/sumA;
	      ave_ccca = sum_ccca/sumA;

	      ave_gs = sum_gs/sumA;
	      ave_gm = sum_gm/sumA;
	    }
	  else
	    {
	      ave_D13 = undef;
	      ave_D13_long = undef;
	      ave_csca = undef;
	      ave_cica = undef;
	      ave_ccca = undef;
	      ave_gs = undef;
	      ave_gm = undef;
	    }

	  // Flux of 13C in micromoles of 13C m-2 s-1

	  // it should relate to Rresp Reco + A (Rair/(1 + D13)

	  fc_13C = sum13C+soil.resp_13;

	  // compute various integrated isotope factors

	  if (can_ps_mol > 0)
	    {
	      for (j=1; j<=jtot; j++)
		{
		  // compute flux weighted discrimination ratio for sun and shaded leaves
		  sun_A = __max(prof.sun_A[j], 0.);
		  shd_A = __max(prof.shd_A[j], 0.);
		  prof.Rresp[j] = (sun_A*solar.prob_beam[j]*prof.Rplant_sun[j] +
					  shd_A*solar.prob_shd[j]*prof.Rplant_shd[j]);
		  prof.Rresp[j] /= (sun_A*solar.prob_beam[j] +
					   shd_A*solar.prob_shd[j]);
		  // sum the isotope ratio to define the respiratory signal based on the carbon that is assimilated
		  if ((shd_A > 0) && (sun_A >0))
		    {
		      prof.Rresp_sum[j] += prof.Rresp[j];
		      prof.cnt_Rresp[j]++;
		    }
		}
	    }
	  else
	    {
	      for (j=1; j<=jtot; j++)
		prof.Rresp[j] = 0.;
	    }
	      
	  // compute isotope flux

	  // computing profile of recycled CO2

	  // theta=(Cf (df-dR)- Ca(da-dR))/(D Cf - (da-dR)Ca)

	  ave_daC13 = 0.; // average del 13C of the air in the canopy
	  //MC ave_Ca_f = 0.; // average CO2 concentration in the forest
	  //MC ave_d13resp = 0.; // ave del 13C of respiring plant

	  for (j=1; j<=jtot; j++)
	    {
	      ave_daC13 += prof.d13Cair[j];
	      //printf("%10.7f ", prof.d13Cair[j]);
	      //MC ave_Ca_f +=prof.co2_air[j];
	      prof.d13Cplant[j] = DELTA1000(prof.Rresp[j], 1., Rpdb_12C);
	      //MC ave_d13resp += prof.d13Cplant[j];
	    }
	  /* printf("\n"); */

	  //MC ave_d13resp /= jtot;
	  ave_daC13 /= jtot;
	  //MC ave_Ca_f /= jtot;

	  // my only concern is weighting plants and soil correctly.

	  // compute ave_d13resp from a simple Keeling plot

	  // Feb 28, 2002 having problems with recyling, getting negative
	  // March 28, 2002 found algebraic error

	  //MC recycle = ave_Ca_f*(ave_daC13-ave_d13resp) - input.co2air*(input.d13CO2-ave_d13resp);
	      
	  //MC recycle /= ave_Ca_f*ave_D13 - input.co2air*(input.d13CO2-ave_d13resp);

	  // Reco = Rsoil + Rbole

	  // Bowling et al. 2001 Global Change Biology 7, 127-145

	  // Fisotope = d13resp Reco +(d13air-D13)A

	  // Need to adjust the signs with the direction of the fluxes to do the accounting correctly
	  // plus up like respiration, negative down like Ps

	  // Ecosystem respiration contains soil and bole respiration.

	  //MC Fisotope = -26.32 * (soil.respiration_mole+bole.respiration_mole) + (ave_daC13-ave_D13)*(-can_ps_mol);
	} // 13C

	  // calculate canopy albedo

      flux.albedo = (solar.par_up[jtot1]/4.6+solar.nir_up[jtot1])
	/ ((solar.beam_flux_par[jtot1] + solar.par_down[jtot1])/4.6
	   + solar.beam_flux_nir[jtot1] + solar.nir_dn[jtot1]);

      if (set_switch.isoprene == 1)
	{
	  // compute isoprene flux density
	  ISOPRENE_CANOPY_FLUX(&isoprene_efflux);
	}

    endday:

      // Daily sums of hourly data

      sumevap += sumle; // latent heat flux
      sumsens += sumh; // sensible heat flux
      sumfc += fc_mol; // CO2 flux, micromoles m-2 s-1
      sumpar += input.parin; // PAR
      sumnet += netrad; // Net radiation
      sumbole += bole.respiration_mole; // bole respiration
      sumsoil += soil.respiration_mole;
      sumresp += canresp; // canopy respiration
      sumta += tleaf_mean; // mean leaf temperature
      sumgs += sumksi; // canopy stomatal conductance
      if (set_switch.d13c == 1)
	sumF13C += fc_13C; // net 13C flux

      daycnt += 1;

      sumps += can_ps_mol; // canopy photosynthesis
      if (can_gpp > 0)
	{
	  sumgpp += can_gpp; // GPP
	  sumTleaf += tleaf_mean; // daytime leaf temperature
	  cntps +=1;
	}
      if (set_switch.d13c == 1)
	{
	  if (can_gpp > 0)
	    {
	      sumD13C += ave_D13*can_ps_mol; // isotopic discrimination
	      sumD13C_long_day += ave_D13_long*can_ps_mol;
	    }
	}

      // write to output structure
      output.netrad = netrad;
      output.sumrn = sumrn;
      output.sumlout = sumlout;
      output.sumh = sumh;
      output.sumle = sumle;
      output.can_ps_mol = can_ps_mol;
      output.can_gpp = can_gpp;
      output.c_transpiration_mole = c_transpiration_mole;
      output.canresp = canresp;
      output.sumksi = sumksi;
      output.tleaf_mean = tleaf_mean;
      output.tavg_sun = tavg_sun;
      output.tavg_shd = tavg_shd;
      output.ave_cica = ave_cica;
      output.ave_gs = ave_gs;
      if (input.parin > 5.)
	output.diff_par = solar.par_diffuse/input.parin;
      else
	output.diff_par = undef;
      output.rnet_soil = rnet_soil;
      for (j=1; j<=jtot; j++)
	output.dLAIdz[j] = prof.dLAIdz[j];
      if (set_switch.d13c == 1)
	{
	  output.ave_D13 = ave_D13;
	  output.ave_D13_long = ave_D13_long;
	  output.ave_daC13 = ave_daC13;
	}
      /*
	WRITE_DUMMY(2,
	sumle-soil.evap,
	temp1
	);
      */

      // DISK OUTPUT

      fprintf(fptr6,"%ld,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,"
	      "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,"
	      "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,"
	      "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,"
	      "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,"
	      "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,"
	      "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,"
	      "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,"
	      "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,"
	      "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,"
	      "%g,%g,%g,%g,%g,%g,%g\n",
	      time_var.daytime, netrad, sumrn, sumh, sumle,
	      can_ps_mol, can_gpp, c_transpiration_mole, canresp,soil.respiration_mole,
	      bole.respiration_mole, soil.gsoil, sumksi, tleaf_mean, tavg_sun,
	      tavg_shd, ave_D13, ave_D13_long, ave_D13_a, ave_D13_ab, ave_D13_asal, ave_D13_b,
	      ave_csca, ave_cica, ave_ccca, ave_gm, ave_gs,
	      ave_daC13, isoprene_efflux, output.diff_par, input.parin,input.ta,
	      input.wnd, met.ustar, met.rhova_g, met.zl, met.press_Pa,
	      met.relative_humidity, input.co2air, prof.tair[40], 
	      prof.co2_air[40], prof.vpd_air[40],
	      prof.tleaf[40], prof.sun_A[40], prof.sun_gs_mol[40], 
	      prof.sun_rbco2[40], prof.sun_cs[40],
	      prof.sun_ci[40], prof.sun_cc[40], prof.sun_wj[40], 
	      prof.sun_wc[40], solar.quantum_sun[40],
	      prof.sun_rh[40], prof.sun_vpd[40],
	      prof.dHdz[40], prof.dLEdz[40][1],
	      prof.dPsdz[40], prof.tleaf[40], prof.D13C[40],
	      prof.D13C_long[40], prof.D13C_a[40], prof.D13C_ab[40], prof.D13C_asal[40],
	      prof.D13C_b[40], prof.cica[40], prof.dStomCondz[40], prof.isopreneflux[40],
	      solar.beam_flux_par[40], solar.par_down[40],
	      solar.beam_flux_nir[40], solar.nir_dn[40],
	      prof.dHdz[30], prof.dLEdz[30][1], prof.dPsdz[30],
	      prof.tleaf[30], prof.D13C[30],
	      prof.D13C_long[30], prof.D13C_a[30], prof.D13C_ab[30], prof.D13C_asal[30],
	      prof.D13C_b[30], prof.cica[30], prof.dStomCondz[30],prof.isopreneflux[30],
	      solar.beam_flux_par[30], solar.par_down[30],solar.beam_flux_nir[30],solar.nir_dn[30],
	      prof.dHdz[20], prof.dLEdz[20][1], prof.dPsdz[20],
	      prof.tleaf[20], prof.D13C[20],
	      prof.D13C_long[20], prof.D13C_a[20], prof.D13C_ab[20], prof.D13C_asal[20],
	      prof.D13C_b[20], prof.cica[20], prof.dStomCondz[20],prof.isopreneflux[20],
	      solar.beam_flux_par[20], solar.par_down[20],solar.beam_flux_nir[20],solar.nir_dn[20],
	      prof.dHdz[10], prof.dLEdz[10][1], prof.dPsdz[10], 
	      prof.tleaf[10], prof.D13C[10],
	      prof.D13C_long[10], prof.D13C_a[10], prof.D13C_ab[10], prof.D13C_asal[10],
	      prof.D13C_b[10], prof.cica[10], prof.dStomCondz[10],prof.isopreneflux[10],
	      solar.beam_flux_par[10], solar.par_down[10],solar.beam_flux_nir[10],solar.nir_dn[10],
	      prof.dHdz[1], prof.dLEdz[1][1], prof.dPsdz[1], prof.tleaf[1], prof.D13C[1],
	      prof.D13C_long[1], prof.D13C_a[1], prof.D13C_ab[1], prof.D13C_asal[1],
	      prof.D13C_b[1], prof.cica[1], prof.dStomCondz[1],prof.isopreneflux[1],
	      solar.beam_flux_par[1], solar.par_down[1],solar.beam_flux_nir[1],
	      solar.nir_dn[1]);

      if (set_switch.d13c == 1)
	fprintf(fptr21, "%ld,%g,%g,%g,%g,%g\n",
		time_var.daytime, output.ave_D13, output.ave_D13_long,
		output.ave_daC13, prof.D13C_long[40],
		prof.D13C[40]);

      fprintf(fptr7, "%ld,%g,%g,%g,%g,%g,%g,"
	      "%g,%g,%g,%g,%g,%g,%g,%g,%g,"
	      "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
	      time_var.daytime, output.rnet_soil, soil.heat,
	      soil.evap, soil.tsfc,
	      soil.T_base, soil.T_15cm, flux.c_transpiration[1],
	      flux.s_evap[1], prof.throughfall[1][1],
	      soil.soil_mm, soil.qinfl[1], soil.qdrai[1],
	      output.c7,soil.qtran[1],flux.surfrun[1],
	      soil.T_soil[1],soil.T_soil[2],soil.T_soil[3],soil.T_soil[4],soil.T_soil[5],
	      soil.T_soil[6],soil.T_soil[7],soil.T_soil[8],soil.T_soil[9],soil.T_soil[10]);

      fprintf(fptr11, "%ld, %g\n",
	      time_var.daytime, soil.soil_mm_50);

      fprintf(fptr13, "%ld,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e,"
	      "%16.10e\n",
	      time_var.daytime,
	      wiso.lost[1],
	      flux.litterevap[1],
	      flux.soilevap[1],
	      flux.s_evap[1],
	      flux.c_evaporation[1],
	      flux.c_transpiration[1],
	      flux.c_evapotranspiration[1],
	      flux.evapotranspiration[1],
	      prof.rhov_air[1][1],
	      prof.rhov_air[jtot3][1],
	      input.ppt[1],
	      prof.throughfall[1][1],
	      flux.surfrun[1],
	      soil.qdrai[1],
	      soil.theta_l[1],
	      soil.theta[1][1],
	      soil.theta[2][1],
	      soil.theta[3][1],
	      soil.theta[4][1],
	      soil.theta[5][1],
	      soil.theta[6][1],
	      soil.theta[7][1],
	      soil.theta[8][1],
	      soil.theta[9][1],
	      soil.theta[10][1]
	      );

      if (set_switch.wiso == 1)
	{
	  fprintf(fptr16, "%ld,"
		  "%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e\n",
		  time_var.daytime,
		  soil.rxylem[2],soil.rxylem[3],soil.rxylem[4],
		  prof.rvapour[1][2],prof.rvapour[1][3],prof.rvapour[1][4]
		  );

	  fprintf(fptr19, "%ld,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e,"
		  "%16.10e,%16.10e,%16.10e,%16.10e\n",
		  time_var.daytime,
		  wiso.lost[1],wiso.lost[2],wiso.lost[3],wiso.lost[4],
		  flux.litterevap[1],flux.litterevap[2],flux.litterevap[3],flux.litterevap[4],
		  flux.soilevap[1],flux.soilevap[2],flux.soilevap[3],flux.soilevap[4],
		  flux.s_evap[1],flux.s_evap[2],flux.s_evap[3],flux.s_evap[4],
		  flux.c_evaporation[1],flux.c_evaporation[2],flux.c_evaporation[3],flux.c_evaporation[4],
		  flux.c_transpiration[1],flux.c_transpiration[2],flux.c_transpiration[3],flux.c_transpiration[4],
		  flux.c_evapotranspiration[1],flux.c_evapotranspiration[2],flux.c_evapotranspiration[3],flux.c_evapotranspiration[4],
		  flux.evapotranspiration[1],flux.evapotranspiration[2],flux.evapotranspiration[3],flux.evapotranspiration[4],
		  prof.rhov_air[1][1], prof.rhov_air[1][2], prof.rhov_air[1][3], prof.rhov_air[1][4],
		  prof.rhov_air[jtot3][1], prof.rhov_air[jtot3][2], prof.rhov_air[jtot3][3], prof.rhov_air[jtot3][4],
		  input.ppt[1], input.ppt[2], input.ppt[3], input.ppt[4],
		  prof.throughfall[1][1],prof.throughfall[1][2],prof.throughfall[1][3],prof.throughfall[1][4],
		  flux.surfrun[1],flux.surfrun[2],flux.surfrun[3],flux.surfrun[4],
		  soil.qdrai[1],soil.qdrai[2],soil.qdrai[3],soil.qdrai[4],
		  soil.theta_l[1],soil.theta_l[2],soil.theta_l[3],soil.theta_l[4],
		  soil.theta[1][1],soil.theta[1][2],soil.theta[1][3],soil.theta[1][4],
		  soil.theta[2][1],soil.theta[2][2],soil.theta[2][3],soil.theta[2][4],
		  soil.theta[3][1],soil.theta[3][2],soil.theta[3][3],soil.theta[3][4],
		  soil.theta[4][1],soil.theta[4][2],soil.theta[4][3],soil.theta[4][4],
		  soil.theta[5][1],soil.theta[5][2],soil.theta[5][3],soil.theta[5][4],
		  soil.theta[6][1],soil.theta[6][2],soil.theta[6][3],soil.theta[6][4],
		  soil.theta[7][1],soil.theta[7][2],soil.theta[7][3],soil.theta[7][4],
		  soil.theta[8][1],soil.theta[8][2],soil.theta[8][3],soil.theta[8][4],
		  soil.theta[9][1],soil.theta[9][2],soil.theta[9][3],soil.theta[9][4],
		  soil.theta[10][1],soil.theta[10][2],soil.theta[10][3],soil.theta[10][4]
		  );
	}

      // output profiles for only specified periods
      if (time_var.daytime >= start_profiles && time_var.daytime <= end_profiles)
	{
	  for (j=1; j<=jtot3; j++)
	    {
	      fprintf(fptr9, "%ld,%i,"
		      "%g,%16.10e,%8.3f\n",
		      time_var.daytime, j, prof.tair[j], prof.rhov_air[j][1],
		      prof.co2_air[j]);
	      if (set_switch.d13c == 1)
		fprintf(fptr22, "%ld,%i,"
			"%12.8f,%10.5f,%8.3f\n",
			time_var.daytime, j, prof.R13_12_air[j],
			prof.d13Cair[j], prof.c13cnc[j]);

	      if (set_switch.wiso == 1)
		fprintf(fptr17, "%ld,%i,"
			"%16.10e,%16.10e,%16.10e\n",
			time_var.daytime, j,
			prof.rhov_air[j][2], prof.rhov_air[j][3], prof.rhov_air[j][4]
			);
	    }
	  for (j=1; j<=jtot; j++)
	    {
	      fprintf(fptr10,"%ld,%i,%g,%g,%g,%g,%g,%g,"
		      "%g,%g,%g,%g,%g,%g,%g,%g,%g,"
		      "%g,%g,%g,%g,%g,%g,%g,%g,%g,"
		      "%g,%g,%g,%g,%g,%g,%g,%g,%g,"
		      "%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
		      time_var.daytime,j,prof.dHdz[j],prof.dLEdz[j][1],
		      prof.dLEdz_sun[j],prof.dLEdz_shd[j], prof.sun_A[j], prof.shd_A[j],
		      prof.dGPPdz[j],prof.dGPPdz_sun[j],prof.dGPPdz_shd[j],
		      prof.dPsdz[j],prof.dPsdz_sun[j],prof.dPsdz_shd[j],
		      prof.dRESPdz[j],prof.dRESPdz_sun[j],prof.dRESPdz_shd[j],
		      solar.quantum_sun[j], solar.quantum_shd[j], prof.tleaf[j],
		      prof.sun_tleaf[j], prof.shd_tleaf[j],
		      solar.prob_beam[j], solar.prob_shd[j], prof.dLAIdz[j], 
		      prof.isopreneflux[j],
		      prof.dStomCondz[j]*1E3,prof.dStomCondz_sun[j]*1E3,prof.dStomCondz_shd[j]*1E3, // conductance in mmol m-2 s-1
		      prof.D13C_long[j], prof.sun_D13_long[j],prof.shd_D13_long[j],
		      prof.D13C[j], prof.sun_D13[j],prof.shd_D13[j],
		      prof.sun_D13_ab[j], prof.sun_D13_a[j], prof.sun_D13_asal[j], 
		      prof.sun_D13_b[j],
		      prof.shd_D13_ab[j], prof.shd_D13_a[j], prof.shd_D13_asal[j], 
		      prof.shd_D13_b[j],
		      input.parin
		      );
		  
	      //if (set_switch.d13c == 1)
	      //  fprintf(fptr23, "%ld,%i\n",
	      //	    time_var.daytime, j);

	      if (set_switch.wiso == 1)
		fprintf(fptr18, "%ld,%i,"
			"%16.10e,%16.10e,%16.10e,"
			"%16.10e,%16.10e,%16.10e,"
			"%16.10e,%16.10e,%16.10e,"
			"%16.10e,%16.10e,%16.10e,"
			"%16.10e,%16.10e,%16.10e,"
			"%16.10e,%16.10e,%16.10e,"
			"%16.10e,%16.10e,%16.10e,"
			"%16.10e,%16.10e,%16.10e,"
			"%16.10e,%16.10e,%16.10e\n",
			time_var.daytime, j,
			prof.dLEdz[j][2], prof.dLEdz[j][3], prof.dLEdz[j][4],
			prof.sun_craig[j][2], prof.sun_craig[j][3], prof.sun_craig[j][4],
			prof.shd_craig[j][2], prof.shd_craig[j][3], prof.shd_craig[j][4],
			prof.sun_leafwater_e[j][2], prof.sun_leafwater_e[j][3], prof.sun_leafwater_e[j][4],
			prof.shd_leafwater_e[j][2], prof.shd_leafwater_e[j][3], prof.shd_leafwater_e[j][4],
			prof.sun_leafwater[j][2], prof.sun_leafwater[j][3], prof.sun_leafwater[j][4],
			prof.shd_leafwater[j][2], prof.shd_leafwater[j][3], prof.shd_leafwater[j][4],
			prof.sun_rtrans[j][2], prof.sun_rtrans[j][3], prof.sun_rtrans[j][4],
			prof.shd_rtrans[j][2], prof.shd_rtrans[j][3], prof.shd_rtrans[j][4]
			);
	    }
	} // output profiles

      /* printf("T: %20.14f %20.14f\n", prof.shd_tleaf_filter[1], prof.shd_tleaf_filter[40]); */

      if (lastin == 1) doit = 0; // interrupt doit loop

    } // end of doit loop

  // Close files and sum up

  fclose(fptr1);   // met in
  fclose(fptr4);   // dispersion matrix
  fclose(fptr6);   // season out
  fclose(fptr7);   // soil out
  fclose(fptr8);   // daily out
  fclose(fptr9);   // profile out
  fclose(fptr10);  // flux profile out
  fclose(fptr11);  // optimise
  fclose(fptr12);  // dummy out
  fclose(fptr13);  // h2o soil out
  if (set_switch.wiso == 1) fclose(fptr14); // wiso in
  if (extra_nate == 1) fclose(fptr15);      // lai in
  if (set_switch.wiso == 1) 
    {
      fclose(fptr16); // wiso h2o leaf out
      fclose(fptr17); // wiso profile out 
      fclose(fptr18); // wiso flux profile out
      fclose(fptr19); // wiso h2o soil out
    }
  if (set_switch.d13c == 1)
    {
      fclose(fptr20); // 13c daily
      fclose(fptr21); // 13c season
      fclose(fptr22); // 13c profile
      //fclose(fptr23); // 13cprofile flux
    }


  end = clock();
  elapsed = ((double)(end-start))/CLOCKS_PER_SEC;
  printf("\n\nStop Canoak.\n");
  printf("\nDone %i time steps. Energy and carbon balance converged %i times = %i%%.\n",
         totalcount, totalcount-count40, (100*(totalcount-count40))/totalcount);
  if (dtotalcount > 0 && ntotalcount > 0)
    printf("\nConverged %i times during daytime = %i%% and %i times during nighttime = %i%%.\n",
           dtotalcount-dcount40, (100*(dtotalcount-dcount40))/dtotalcount,
           ntotalcount-ncount40, (100*(ntotalcount-ncount40))/ntotalcount);
  printf("\nSeconds elapsed since start: %f.\n", elapsed);
  printf("Output written to %s\n", outputfolder);

  return (!ok);
} // end of main function

// ----------------------------------------------------
// Subroutines and Functions
// ----------------------------------------------------

// ----------------------------------------------------
void SKIP_INPUT()
{
  fpos_t pos;
  long daytime;
  int dayy, flag;
  float hhrr, ta, rglobal, parin, pardif, ea, wnd, ppt;
  float co2air, press_mb, tsoil, soilmoisture, d13CO2, d18CO2;

  fpos_t pos1, pos2;
  int dayy1, dayy2;
  float hhrr1, d18o, d2h, d18o2;
  float hhrr2, lai1, lai2;

  fgetpos(fptr1, &pos);
  fscanf(fptr1,"%i %g %g %g %g %g %g %g %g %g %g %g %g %i %g %g\n"
         , &dayy, &hhrr, &ta, &rglobal, &parin, &pardif
         , &ea, &wnd, &ppt, &co2air, &press_mb, &tsoil
         , &soilmoisture, &flag, &d13CO2, &d18CO2);
  daytime = dayy*10000 + (int) (hhrr*100.); // define daytime

  // water isotopes
  if (set_switch.wiso == 1)
    {
      fgetpos(fptr14, &pos1);
      fscanf(fptr14,"%i %g %g %g %g\n"
	     , &dayy1, &hhrr1, &d18o, &d2h, &d18o2);
    }
  if (extra_nate == 1)
    {
      // lai
      fgetpos(fptr15, &pos2);
      fscanf(fptr15,"%i %g %g %g\n"
	     , &dayy2, &hhrr2, &lai1, &lai2);
    }

  while (daytime < start_run && !feof(fptr1))
    {
     
      fgetpos(fptr1, &pos);
      fscanf(fptr1,"%i %g %g %g %g %g %g %g %g %g %g %g %g %i %g %g\n"
             , &dayy, &hhrr, &ta, &rglobal, &parin, &pardif
             , &ea, &wnd, &ppt, &co2air, &press_mb, &tsoil
             , &soilmoisture, &flag, &d13CO2, &d18CO2);
      daytime = dayy*10000 + (int) (hhrr*100.); // define daytime
      // water isotopes
      if (set_switch.wiso == 1)
	{
	  if (feof(fptr14)) puts("Problem with h2o iso input in skip_input.");
	  fgetpos(fptr14, &pos1);
	  fscanf(fptr14,"%i %g %g %g %g\n"
		 , &dayy1, &hhrr1, &d18o, &d2h, &d18o2);
	}
      if (extra_nate == 1)
	{
	  // lai
	  if (feof(fptr15)) puts("Problem with extra lai input in skip_input.");
	  fgetpos(fptr15, &pos2);
	  fscanf(fptr15,"%i %g %g %g\n"
		 , &dayy2, &hhrr2, &lai1, &lai2);
	}
    }

  fsetpos(fptr1, &pos);
  if (set_switch.wiso == 1) fsetpos(fptr14, &pos1);
  if (extra_nate == 1) fsetpos(fptr15, &pos2);

  return;
}


// ----------------------------------------------------
void INPUT_DATA()
{
  // input data and check for bad data
  // note that the data were produced in single precision (float)
  // so I had to read them as single precision, otherwise I ingested
  // garbage
  double est=0.;
  int j=0, mc=0;
  //MC20190615 float tmp;
  double tmp;
  long int dt;

  int dayy1=0, dayy2=0;
  //MC20190615 float hhrr1=0., d18o=0., d2h=0., d18ov=0., hhrr2=0.;
  double hhrr1=0., d18o=0., d2h=0., d18ov=0., hhrr2=0.;

  dt = ((long int) time_var.time_step) * 100 / 3600;
  // use input.dayy instead of time_var.days because of several years at once
  time_var.jdold = input.dayy; // identify previous day
      
  //MC20190615 fscanf(fptr1,"%i %g %g %g %g %g %g %g %g %g %g %g %g %ld %g %g\n"
  fscanf(fptr1,"%i %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %ld %lg %lg\n"
	 , &input.dayy, &input.hhrr, &input.ta
	 , &input.rglobal, &input.parin, &input.pardif
	 , &input.ea, &input.wnd, &tmp
	 , &input.co2air, &input.press_mb, &input.tsoil
	 , &input.soilmoisture, &input.flag, &input.d13CO2
	 , &input.d18CO2);
  //MC20190615 input.ppt[1] = (double) tmp;
  input.ppt[1] = tmp;
  /* printf("RI01.01 %d %20.14f %20.14f\n", input.dayy, input.hhrr, input.ta); */
  /* printf("RI01.02 %20.14f %20.14f %20.14f\n", input.rglobal, input.parin, input.pardif); */
  /* printf("RI01.03 %20.14f %20.14f %20.14f\n", input.ea, input.wnd, input.ppt[1]); */
  /* printf("RI01.04 %20.14f %20.14f %20.14f\n", input.co2air, input.press_mb, input.tsoil); */
  /* printf("RI01.05 %20.14f %ld %20.14f\n", input.soilmoisture, input.flag, input.d13CO2); */
  /* printf("RI01.06 %20.14f\n", input.d18CO2); */

  // for Nate McDowell's juniper site read LAI instead of diffuse PAR
  if (extra_nate == 1)
    input.lai = input.pardif; // redundant with fptr15

  // preliminary approach to calculate consecutive year
  time_var.year = time_var.year0 + (input.dayy-1)/365;

  time_var.daytime = input.dayy*10000 + (int) (input.hhrr*100); // define daytime
  
  time_var.doy = input.dayy % 365;    
  if (time_var.doy == 0)
    time_var.doy = 365;
  
  time_var.local_time = input.hhrr;

  time_var.days = time_var.doy;

  // compute derived quantities for the model

  met.T_Kelvin = input.ta + TN0; // compute absolute air temperature

  met.rhova_g = input.ea * 2165./met.T_Kelvin; // compute absolute humidity, g m-3

  // limit humidity
  if (met.rhova_g < 0.)
    met.rhova_g = 0.;
  if (met.rhova_g > 30.)
    met.rhova_g = 30.;

  met.rhova_kg = met.rhova_g / 1000.; /* absolute humidity, kg m-3 */

  met.press_kpa = input.press_mb/10.; // air pressure, kPa
  met.press_bars = input.press_mb/1000.; // air pressure, bars
  met.press_Pa = met.press_kpa*1000.; // pressure, Pa

  est = ES(met.T_Kelvin);
  met.relative_humidity = input.ea*10./est; // relative humidity

  // combining gas law constants
  met.pstat273 = rugc / (100000. * met.press_bars);

  for (j=1; j<=jtot; j++)
    {
      // cuticular conductance adjusted for pressure and T, mol m-2 s-1
      sfc_res.gcut = bprime[j] * met.T_Kelvin * met.pstat273;
      // cuticular resistance
      sfc_res.rcuticle[j] = 1.0 / sfc_res.gcut;
    }

  if (fabs(input.co2air) >= 998.) // check for bad CO2 data
    input.co2air = 370.;
  
  /* printf("RI01.07 %20.14f\n", input.parin); */
  if (input.parin < 0.)
    input.parin = 0.; // check for bad par data

  /* printf("RI01.08 %20.14f\n", input.parin); */
  //MC20190625 solar.* set later in code in diffuse_direct_radiation()
  /* if (input.parin <= 0.) // check for night */
  /*   { */
  /*     solar.par_beam = 0.; */
  /*     solar.par_diffuse = 0.; */
  /*     solar.nir_beam = 0.; */
  /*     solar.nir_diffuse = 0.; */
  /*   } */
  /* printf("RI01.09 %20.14f %20.14f %20.14f\n", solar.par_diffuse, solar.nir_beam, solar.nir_diffuse); */
  /* printf("RI01.10 %20.14f\n", solar.ratrad); */

  // set some limits on bad input data to minimize the model from blowing up
  if (solar.ratrad > 0.9 || solar.ratrad < 0.2)
    solar.ratrad=0.5;
  /* printf("RI01.11 %20.14f\n", solar.ratrad); */

  // limit wind speed
  if (input.wnd < 1.)
    input.wnd = 1.;
  if (input.wnd > 10.)
    input.wnd = 10.;

  met.ustar_filter = vonKarman/log((zm-zd)/z0)*input.wnd;

  // limit air temperature
  if (input.ta < -30.)
    input.ta = -30.;
  if (input.ta > 60.)
    input.ta = 60.;

  // air density, kg m-3
  met.air_density = met.press_kpa * mass_air / (rugc * met.T_Kelvin);

  // air density, mole m-3
  met.air_density_mole = met.press_kpa/ (rugc * met.T_Kelvin) * 1000.;

  // water isotopes
  if (set_switch.wiso == 1)
    {
      if (feof(fptr14)) puts("Problem with h2o iso input in input_data.");
      //MC20190615 fscanf(fptr14,"%i %g %g %g %g\n"
      fscanf(fptr14,"%i %lg %lg %lg %lg\n"
	     , &dayy1, &hhrr1, &d18o, &d2h, &d18ov);
      input.dppt[1] = 0.;
      //MC20190615 input.dppt[2] = (double) d18o;
      input.dppt[2] = d18o;
      //MC20190615 input.dppt[3] = (double) d2h;
      input.dppt[3] = d2h;
      input.dppt[4] = 0.;
      //MC20190615 input.dvapour[2] = (double) d18ov;
      input.dvapour[2] = d18ov;
      
      // Test Meteoric water line for deuterium
      input.dppt[3] = 8. * input.dppt[2] + 10.;

      if (wiso.nofracin == 1)
	for (mc=2; mc<=wiso.nwater+1; mc++)
	  input.dppt[mc] = wiso.dtheta[1][mc];

      for (mc=2; mc<=wiso.nwater+1; mc++)
	input.ppt[mc] = input.ppt[1] * INVDELTA1000_H2O(input.dppt[mc], mc);
    }

  if (extra_nate == 1)
    {
      // LAI file
      if (feof(fptr15)) puts("Problem with extra lai input in input_data.");
      //MC20190615 fscanf(fptr15,"%i %g %g %g\n"
      fscanf(fptr15,"%i %lg %lg %lg\n"
	     , &dayy2, &hhrr2, &input.lai_up, &input.lai_down);
      // Test: set LAI to fix value of 1 + 1 m2 m-2
      /* input.lai_up = input.lai_down = 1; */
      input.lai_up   = 1;
      input.lai_down = 1;
      input.lai = input.lai_up + input.lai_down;
      lai = 2;
    }

  if ((time_var.daytime+dt) > end_run)
    lastin = 1;
  
  return;
}


// ----------------------------------------------------
void SOIL_RESPIRATION()
{
  // Computes soil respiration

  double fautleaf = 0.40;
  double ccost = 0.1;
  
  double zass = 0.;   // GPP of entire canopy layer
  double zrd = 0.;   // dark respiration of leaves 
  double zauto = 0.;  // autotrophic respiration of entire canopy
  double ztmp2d = 0.; // maintenace respiration
  int j;

  /* OLD
     soil.base_respiration=0.71; // new value Recosystem=1.34
     soil.respiration_mole = soil.base_respiration * exp((51609. / 8.314) * ((1. / 273.) - 1. / (soil.T_15cm + TN0)));

     // soil wetness factor from the Hanson model, assuming constant and wet soils
     soil.respiration_mole *= 0.86;
  */

  if (extra_nate == 1)
    {
      // Equation of Nate McDowell for Juniper site
      // Take 5cm instead of Nate's 2cm because Canoak's 2cm is too variable
      soil.respiration_mole = 0.096*soil.theta[4][1]*100. + 0.5089; 
      /* printf("R: %20.14f %20.14f\n", soil.respiration_mole, soil.theta[4][1]); */
    }
  else
    {
      // calculate soil respiration based on chamber measurements by Astrd Soe for Hainich site (Soe 2003)
      // soil respiraton = f(Tsoil in 5 cm)

      if (set_switch.soil_resp_temp == 1)
	soil.SR_ref_temp = input.tsoil; // set soil reference temperature to input soil temperature
      else
	soil.SR_ref_temp = soil.T_soil[4]; // set soil reference temperature to modeled soil temperature in 5 cm, depending on switch set in parameter file
  
      // canisotope v3.1
      //soil.respiration_mole = 0.71*exp(0.09*soil.SR_ref_temp);	//factor changed, orignal 0.11
      // from Hainich soil respiration measurements
      soil.respiration_mole = 0.83*exp(0.11*soil.SR_ref_temp);
    }

  
  if (set_switch.bethy_resp == 1)
    {
      /* 
	 Compute soil respiration as autotrophic and heterotrophic respiration based on BETHY (Knorr PhD thesis 1997)
	 Autotrophic is maintenance and growth respiration
	 Rleaf = 0.4 Rmaintenance => Rmaintenance = Rleaf/0.4
	 Rgrowth = 0.25 * NPP = 0.25*(GPP - Rmaintenence - Rgrowth) => Rgrowth = 0.25/(1+0.25) * (GPP - Rmaintenance)
      */
      for (j=1; j<=jtot; j++)
	{
	  zass += prof.dPsdz[j] + prof.dRESPdz[j];
	  zrd += prof.dRESPdz[j];
	}

      ztmp2d = zrd/fautleaf;  // maintenance respiration of entire plants
      zauto = ztmp2d + __max((zass-ztmp2d)*(ccost/(1.+ccost)), 0.); //maintenance + growth respiration

      // soil autotroph respiration = total autotroph - leaf - bole
      soil.respiration_auto = zauto-zrd-bole.respiration_mole;
      soil.respiration_hetero = soil.respiration_mole - soil.respiration_auto;
    }
  else
    {
      // assume autotroph and heterotroph are both 50% of soil respiration
      soil.respiration_auto = soil.respiration_mole/2.;
      soil.respiration_hetero = soil.respiration_mole/2.;
    }

  // convert soilresp to mg m-2 s-1 from umol m-2 s-1
  soil.respiration_mg = soil.respiration_mole * mass_CO2/1000.;

  return;
}


// ----------------------------------------------------
void BOLE_RESPIRATION()
{
  /*
    Compute bole respiration

    Using data of Edwards and Hanson from WBW to
    estimate bole respiration on ground area
    SAI is 1.5. Rbole is the sum of growth and maintenance
    respiration. Growth respiration is
    about 30% of maintenance respiration. Q10 was about 2.4
    for oak and 1.7 for maple. Geometric mean is 2.02.
  */

  double tempbole, air_ratio, stemlay;

  int j;

  tempbole = 0.5*soil.T_15cm + 0.5*input.ta;

  //MC20190625 bole.calc = (tempbole-10.)/(8.314*283*(tempbole+TN0));
  bole.calc = (tempbole-10.)/(rugc*283.*(tempbole+TN0));

  if (time_var.days >= time_var.leafout && time_var.days <= time_var.leaffallcomplete)
    bole.factor = 0.43;
  else
    bole.factor = 0.259;

  // for Nate McDowell's Juniper site
  if (extra_nate == 1)
    bole.factor *= 0.5;

  // bole respiration flux density, micromol m-2 s-1

  bole.respiration_mole = bole.factor * exp(eabole*bole.calc);

  // convert bole.respiration_mole to mg m-2 s-1

  bole.respiration_mg = bole.respiration_mole * mass_CO2/1000.;

  /*
    Divide the bole resp into layers in the stem space and
    subtract from layer ps, mg CO2 m-3 s-1
  */
  air_ratio = mass_air/(met.air_density*mass_CO2);

  for (j=1; j<=jtot; j++)
    {
      stemlay = prof.dPAIdz[j] / pai;
      // bole.layer[j] = bole.respiration_mg * stemlay;

      bole.layer[j] = bole.respiration_mole * stemlay;
      prof.source_co2[j] += bole.layer[j];

      // prof.source_co2[j] *= air_ratio;
    }

  return;
}


// ----------------------------------------------------
void FILENAMES()
{
  // Constructs Filenames and Opens files

  char *dailyfile = "daily_ave";
  char *daily13cfile = "daily_ave_13c";

  char *hourlyfile = "season";
  char *hourly13cfile = "season_13c";

  char *optimisefile = "optimise";

  char *profilefile = "profile_air";
  char *profile13cfile = "profile_air_13c";
  char *profileisofile = "profile_air_wiso";

  char *fluxprofilefile = "profile_fluxes";
  // char *fluxprofile13cfile = "profile_fluxes_13c";
  char *fluxprofileisofile = "profile_fluxes_wiso";

  char *soilfile = "soil";

  char *dummyfile = "dummy";

  char *h2osoilfile = "h2osoil";
  char *h2osoilisofile = "h2osoil_wiso";

  char *h2oleafisofile = "h2oleaf_wiso";

  char outbuff[128], optimisebuff[128], inbuff[128], soilbuff[128], avebuff[128];
  char yrbuff[5], fluxbuff[128], profilebuff[128], dummybuff[128], h2osoilbuff[128], dispbuff[128];
  char isoinbuff[128], laiinbuff[128], h2oleafisobuff[128], h2osoilisobuff[128], fluxisobuff[128], profileisobuff[128];
  char out13cbuff[128], ave13cbuff[128], profile13cbuff[128];//, flux13cbuff[128];
  char minus[2];
  int ok;

  // define variables for reading current time and date
  char current_timedate[128];
  time_t ltime;
  struct tm *today;
  char destination[128];
  //char outputfolder[128]; // define global to print out in main

  puts("Output files:");
  sprintf(minus, "-");
  // read current time and date and save in tmpbuf in yyyymmdd_hhmm format
  time(&ltime);
  today = localtime(&ltime);
  strftime(current_timedate, 128, "%Y%m%d_%H%M", today);

  sprintf(yrbuff, "%i", time_var.year);

  // Open Dispersion Matrix

  sprintf(dispbuff, "%s", dispfile);
  fptr4 = fopen(dispbuff, "r");
  if (fptr4 == NULL) puts("Cannot open dispersion file.\n");

  // create 'model_output' directory if it does not exist
  sprintf(outputfolder, "%s%s%s", disksubdir, "model_output", "/");
#ifdef WIN32
  ok = _mkdir(outputfolder);
#else
  ok = mkdir(outputfolder, S_IRWXU+S_IRWXG+S_IRWXO);
#endif

  // create new folder with current time and date for model output
  if (set_switch.mode == 0)
    sprintf(outputfolder, "%s%s%s%s%s", disksubdir, "model_output", "/", current_timedate, "/");

  // create optimise folder
  if (set_switch.mode == 1)
    sprintf(outputfolder, "%s%s%s%s%s", disksubdir, "model_output", "/", "optimise", "/");

#ifdef WIN32
  ok = _mkdir(outputfolder);
#else
  ok = mkdir(outputfolder, S_IRWXU+S_IRWXG+S_IRWXO);
#endif

  // file name hourly output

  //sprintf(outbuff, "%s%s%s%s%s", outputfolder, hourlyfile, minus, yrbuff, filesuffix);
  sprintf(outbuff, "%s%s%s%s%s", outputfolder, hourlyfile, minus, yrbuff, filesuffix);
  puts(outbuff);

  fptr6=fopen(outbuff, "w");
  if (fptr6 == NULL) puts("Cannot open seasonYR.dat");

  fprintf(fptr6,"daytime,netrn,sumrn,sumh,sumle,"
	  "canps,gpp,transpiration,canresp,soilresp,"
	  "boleresp,gsoil,ksi,Tleaf,Tlsun,"
	  "Tlshade,D13,D13_long,D13_a,D13_ab,D13_asal,D13_b,"
	  "CsCa,CiCa,CcCa,gm,gs,"
	  "ave_dair_13C,isoprene.flux,diff_par,par,ta,u,"
	  "ustar,rhov,zL,press.Pa,relative.humidity,"
	  "CO2air,Tair40,CO2air40,Vpd40,Tleaf40,"
	  "Sun_A40,Sun_gs_mol40,Sun_rbCO240,Sun_cs40,Sun_ci40,"
	  "Sun_cc40,Sun_wj40,Sun_wc40,Quantum_sun40,Sun_rh40,Sun_vpd40,"
	  "H_40,Trans_40,canps_40,Tleaf_40,D13_40,D13C_long_40,D13C_a_40,"
	  "D13C_ab_40,D13C_asl_40,D13C_b_40,cica_40,gs_40,Fisoprene_40,"
	  "PARdirect_40,PARdiffuse_40,NIRdirect_40,NIRdiffuse_40,"
	  "H_30,Trans_30,canps_30,Tleaf_30,D13_30,D13C_long_30,D13C_a_30,"
	  "D13C_ab_30,D13C_asl_30,D13C_b_30,cica_30,gs_30,Fisoprene_30,"
	  "PARdirect_30,PARdiffuse_30,NIRdirect_30,NIRdiffuse_30,"
	  "H_20,Trans_20,canps_20,Tleaf_20,D13_20,D13C_long_20,D13C_a_20,"
	  "D13C_ab_20,D13C_asl_20,D13C_b_20,cica_20,gs_20,Fisoprene_20,"
	  "PARdirect_20,PARdiffuse_20,NIRdirect_20,NIRdiffuse_20,"
	  "H_10,Trans_10,canps_10,Tleaf_10,D13_10,D13C_long_10,D13C_a_10,"
	  "D13C_ab_10,D13C_asl_10,D13C_b_10,cica_10,gs_10,Fisoprene_10,"
	  "PARdirect_10,PARdiffuse_10,NIRdirect_10,NIRdiffuse_10,"
	  "H_1,Trans_1,canps_1,Tleaf_1,D13_1,D13C_long_1,D13C_a_1,D13C_ab_1,"
	  "D13C_asl_1,D13C_b_1,cica_1,gs_1,Fisoprene_1,PARdirect_1,"
	  "PARdiffuse_1,NIRdirect_1,NIRdiffuse_1\n");

  // file name optimse output

  sprintf(optimisebuff, "%s%s%s%s%s", outputfolder, optimisefile, minus, yrbuff, filesuffix);
  puts(optimisebuff);

  fptr11 = fopen(optimisebuff, "w");
  if (fptr11 == NULL) puts("Cannot open optimseYR.dat");

  fprintf(fptr11, "daytime,soil_mm_50\n");

  // file name soil

  sprintf(soilbuff, "%s%s%s%s%s", outputfolder, soilfile, minus, yrbuff, filesuffix);
  puts(soilbuff);

  fptr7 = fopen(soilbuff, "w");
  if (fptr7 == NULL) puts("Cannot open soilYR.dat.\n");

  fprintf(fptr7, "daytime,netrn,soilh,soille,soil.tsfc,"
          "soil.T.base,soil.T.15cm,flux.transp,flux.evap,prof.throughfall,"
          "soil.soil_mm,soil.qinfl,soil.qdrai,soil.gsoil,soil.qtran,soil.surfrun,"
          "Tsoil.1,Tsoil.2,Tsoil.3,Tsoil.4,Tsoil.5,Tsoil.6,Tsoil.7,Tsoil.8,Tsoil.9,Tsoil.10\n");

  // Name file for daily average output

  sprintf(avebuff, "%s%s%s%s%s", outputfolder, dailyfile, minus, yrbuff, filesuffix);
  puts(avebuff);

  fptr8=fopen(avebuff,"w");
  if (fptr8 == NULL) puts("Cannot open aveYR.dat.\n");

  fprintf(fptr8, "Day,Avg_FC,Avg_EVAP,AVG_H,Avg_PAR,Avg_RNET,lai,Avg_PS,"
          "Ave_Resp,Avg_BOLE,Avg_SOIL,Avg_TLeaf,Avg_Gs,Tleaf_day\n");

  // profiles

  // "profile_air"

  sprintf(profilebuff, "%s%s%s%s%s", outputfolder, profilefile, minus, yrbuff, filesuffix);
  puts(profilebuff);

  fptr9 = fopen(profilebuff, "w");
  if (fptr9 == NULL) puts("Cannot open profileYR.dat.\n");

  fprintf(fptr9, "Time,i,tair,qair,co2\n");

  // "profile_fluxes"

  sprintf(fluxbuff, "%s%s%s%s%s", outputfolder, fluxprofilefile, minus, yrbuff, filesuffix);
  puts(fluxbuff);

  fptr10 = fopen(fluxbuff, "w");
  if (fptr10 == NULL) puts("Cannot open fluxprofileYR.dat.\n");

  fprintf(fptr10,"daytime,i,dHdz,dLEdz,dLEdz.sun,dLEdz.shd,sun_A,shd_A,"
	  "dGPPdz,dGPPdz.sun,dGPPdz.shd,dPsdz,dPsdz.sun,dPsdz.shd,dRESPdz,"
	  "dRESPdz.sun,dRESPdz.shd,"
	  "PARdirect,PARdiffuse,Tleaf,Tleaf_sun,Tleaf_shd,"
	  "Beam,Nonbeam,LAI,Isopreneflux,"
	  "gs,gs.sun,gs.shd,"
	  "D13Clong,D13Clong.sun,D13Clong.shd,D13C,D13C.sun,D13C.shd,"
	  "D13Ca.sun,D13Cab.sun,D13Casal.sun,D13Cb.sun,"
	  "D13Ca.shd,D13Cab.shd,D13Casal.shd,D13Cb.shd,"
	  "PARin\n");

  // Dummy file with 50 columns

  sprintf(dummybuff, "%s%s%s%s%s", outputfolder, dummyfile, minus, yrbuff, filesuffix);
  puts(dummybuff);

  fptr12 = fopen(dummybuff, "w");
  if (fptr12 == NULL) puts("Cannot open dummy file");

  fprintf(fptr12, "Time,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,"
          "dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,"
          "dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,"
          "dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum\n");

  //  Soil water file
  sprintf(h2osoilbuff, "%s%s%s%s%s", outputfolder, h2osoilfile, minus, yrbuff, filesuffix);
  puts(h2osoilbuff);

  fptr13 = fopen(h2osoilbuff, "w");
  if (fptr13 == NULL) puts("Cannot open h2osoil file");

  fprintf(fptr13, "Time,"
              "wiso.lost1,"
              "flux.litterevap1,"
              "flux.soilevap1,"
              "flux.evap1,"
              "flux.c_evaporation1,"
              "flux.c_transpiration1,"
              "flux.c_evapotranspiration1,"
              "flux.evapotranspiration1,"
              "prof.rhov_air11,"
              "prof.rhov_air1201,"
              "input.ppt1,"
              "prof.throughfall1,"
              "flux.surfrun1,"
              "soil.qdrai1,"
              "soil.thetal1,"
              "soil.theta01,"
	      "soil.theta11,"
	      "soil.theta21,"
	      "soil.theta31,"
	      "soil.theta41,"
	      "soil.theta51,"
	      "soil.theta61,"
	      "soil.theta71,"
	      "soil.theta81,"
	      "soil.theta91\n"
              );

  // 13CO2 isotope files

  if (set_switch.d13c == 1)
    {
      // 13C daily average output
      sprintf(ave13cbuff, "%s%s%s%s%s", outputfolder, daily13cfile, minus, yrbuff, filesuffix);
      puts(ave13cbuff);

      fptr20 = fopen(ave13cbuff,"w");
      if (fptr20 == NULL) puts("Cannot open 13C aveYR.dat.\n");

      fprintf(fptr20, "Day,Avg_Fc_13C,Avg_D13_day\n");

      // 13C hourly output: season
      sprintf(out13cbuff, "%s%s%s%s%s", outputfolder, hourly13cfile, minus, yrbuff, filesuffix);
      puts(out13cbuff);

      fptr21 = fopen(out13cbuff, "w");
      if (fptr21 == NULL) puts("Cannot open 13C seasonYR.dat");

      fprintf(fptr21, "daytime,D13,D13_long,"
          "ave_dair_13C,D13Clong40,D13C40\n");

      // profiles

      // 13C profile air

      sprintf(profile13cbuff, "%s%s%s%s%s", outputfolder, profile13cfile, minus, yrbuff, filesuffix);
      puts(profile13cbuff);

      fptr22 = fopen(profile13cbuff, "w");
      if (fptr22 == NULL) puts("Cannot open 13C profileYR.dat.\n");

      fprintf(fptr22, "Time,i,R13C,del13C,13C\n");

      // 13C profile fluxes

      //sprintf(flux13cbuff, "%s%s%s%s%s", outputfolder, flux13cprofilefile, minus, yrbuff, filesuffix);
      //puts(flux13cbuff);

      //fptr23 = fopen(flux13cbuff, "w");
      //if (fptr23 == NULL) puts("Cannot open 13C fluxprofileYR.dat.\n");
      
      //fprintf(fptr23, "daytime, i\n");
    } // 13C

  // Water isotope files

  if (set_switch.wiso == 1)
    {
      // Iso leaf water
      sprintf(h2oleafisobuff, "%s%s%s%s%s", outputfolder, h2oleafisofile, minus, yrbuff, filesuffix);
      puts(h2oleafisobuff);

      fptr16 = fopen(h2oleafisobuff,"w");
      if (fptr16 == NULL) puts("Cannot open h2oleaf_wiso file");

      fprintf(fptr16, "Time,"
              "soil.rxylem2,soil.rxylem3,soil.rxylem4,"
              "prof.rvapour12,prof.rvapour13,prof.rvapour14\n"
              );

      // Iso profile
      sprintf(profileisobuff, "%s%s%s%s%s", outputfolder, profileisofile, minus, yrbuff, filesuffix);
      puts(profileisobuff);

      fptr17 = fopen(profileisobuff, "w");
      if (fptr17 == NULL) puts("Cannot open profile_wisoYR.dat.\n");
      fprintf(fptr17, "Time,i,"
	      "qair2,qair3,qair4\n"
	      );

      // Iso flux profile
      sprintf(fluxisobuff, "%s%s%s%s%s", outputfolder, fluxprofileisofile, minus, yrbuff, filesuffix);
      puts(fluxisobuff);

      fptr18 = fopen(fluxisobuff, "w");
      if (fptr18 == NULL) puts("Cannot open fluxprofile_wisoYR.dat.\n");

      fprintf(fptr18, "daytime,i,"
	      "dLEdz2,dLEdz3,dLEdz4,"
	      "sun_craig2,sun_craig3,sun_craig4,"
	      "shd_craig2,shd_craig3,shd_craig4,"
	      "sun_leafwater_e2,sun_leafwater_e3,sun_leafwater_e4,"
	      "shd_leafwater_e2,shd_leafwater_e3,shd_leafwater_e4,"
	      "sun_leafwater2,sun_leafwater3,sun_leafwater4,"
	      "shd_leafwater2,shd_leafwater3,shd_leafwater4,"
	      "sun_trans2,sun_trans3,sun_trans4,"
	      "shd_trans2,shd_trans3,shd_trans4\n");

      // Iso soil water
      sprintf(h2osoilisobuff, "%s%s%s%s%s", outputfolder, h2osoilisofile, minus, yrbuff, filesuffix);
      puts(h2osoilisobuff);

      fptr19 = fopen(h2osoilisobuff,"w");
      if (fptr19 == NULL) puts("Cannot open h2osoil_wiso file");

      fprintf(fptr19, "Time,"
              "wiso.lost1,wiso.lost2,wiso.lost3,wiso.lost4,"
              "flux.litterevap1,flux.litterevap2,flux.litterevap3,flux.litterevap4,"
              "flux.soilevap1,flux.soilevap2,flux.soilevap3,flux.soilevap4,"
              "flux.s_evap1,flux.s_evap2,flux.s_evap3,flux.s_evap4,"
              "flux.c_evaporation1,flux.c_evaporation2,flux.c_evaporation3,flux.c_evaporation4,"
              "flux.c_transpiration1,flux.c_transpiration2,flux.c_transpiration3,flux.c_transpiration4,"
              "flux.c_evapotranspiration1,flux.c_evapotranspiration2,flux.c_evapotranspiration3,flux.c_evapotranspiration4,"
              "flux.evapotranspiration1,flux.evapotranspiration2,flux.evapotranspiration3,flux.evapotranspiration4,"
              "prof.rhov_air11,prof.rhov_air12,prof.rhov_air13,prof.rhov_air14,"
              "prof.rhov_air1201,prof.rhov_air1202,prof.rhov_air1203,prof.rhov_air1204,"
              "input.ppt1,input.ppt2,input.ppt3,input.ppt4,"
              "prof.throughfall1,prof.throughfall2,prof.throughfall3,prof.throughfall4,"
              "flux.surfrun1,flux.surfrun2,flux.surfrun3,flux.surfrun4,"
              "soil.qdrai1,soil.qdrai2,soil.qdrai3,soil.qdrai4,"
              "soil.thetal1,soil.thetal2,soil.thetal3,soil.thetal4,"
              "soil.theta01,soil.theta02,soil.theta03,soil.theta04,"
	      "soil.theta11,soil.theta12,soil.theta13,soil.theta14,"
	      "soil.theta21,soil.theta22,soil.theta23,soil.theta24,"
	      "soil.theta31,soil.theta32,soil.theta33,soil.theta34,"
	      "soil.theta41,soil.theta42,soil.theta43,soil.theta44,"
	      "soil.theta51,soil.theta52,soil.theta53,soil.theta54,"
	      "soil.theta61,soil.theta62,soil.theta63,soil.theta64,"
	      "soil.theta71,soil.theta72,soil.theta73,soil.theta74,"
	      "soil.theta81,soil.theta82,soil.theta83,soil.theta84,"
	      "soil.theta91,soil.theta92,soil.theta93,soil.theta94\n"
              );
    }

  // name the input file
  sprintf(inbuff, "%s%s%s", input_file_subdir, yrbuff, input_file_suffix);

  fptr1 = fopen(inbuff, "r");
  if (fptr1 == NULL) puts("Cannot open input file.\n");
  //  lines_in_file = 0;
  //  while (fgets(test, 256, fptr1) != NULL) lines_in_file++;
  //  lines_in_file--; // substract header
  //  rewind(fptr1);

  // water isotopes
  if (set_switch.wiso == 1)
    {
      // name the isotope input file
      sprintf(isoinbuff, "%s%swiso%s", input_file_subdir, yrbuff, input_file_suffix);

      fptr14 = fopen(isoinbuff, "r");
      if (fptr14 == NULL) puts("Cannot open wiso input file.\n");
    }

  if (extra_nate == 1)
    {
      // LAI
      // name the lai input file
      sprintf(laiinbuff, "%s%slai%s", input_file_subdir, yrbuff, input_file_suffix);

      fptr15 = fopen(laiinbuff, "r");
      if (fptr15 == NULL) puts("Cannot open lai input file.\n");
    }

  // copy this file, parameter file and input file to output folder
  sprintf(destination, "%s%s", outputfolder, parameter_file);
  FileCopy(parameter_file, destination);

  sprintf(destination, "%s%s", outputfolder, code_file);
  FileCopy(code_file, destination);

  sprintf(destination, "%s%s%s", outputfolder, yrbuff, input_file_suffix);
  FileCopy(inbuff, destination);

  // water isotopes
  if (set_switch.wiso == 1)
    {
      sprintf(destination, "%s%swiso%s", outputfolder, yrbuff, input_file_suffix);
      FileCopy(isoinbuff, destination);
    }

  if (extra_nate == 1)
    {
      // LAI
      sprintf(destination, "%s%slai%s", outputfolder, yrbuff, input_file_suffix);
      FileCopy(laiinbuff, destination);
    }

  puts("");

  return;
}


// ----------------------------------------------------
void G_FUNC_DIFFUSE(){

  /*
    ------------------------------------------------------------
    This subroutine computes the G Function according to
    the algorithms of Lemeur (1973, Agric. Meteorol. 12: 229-247).

    The original code was in Fortran and was converted to C

    This program computs G for each sky sector, as
    needed to compute the transmission of diffuse light

    G varies with height to account for vertical variations
    in leaf angles

    -------------------------------------------------------
  */


  int j, IN, K, KK, I, II, IN1;

  double aden[18], TT[18], PGF[18], sin_TT[18], del_TT[18], del_sin[18];
  double PPP;
  double ang,dang, aang;
  double cos_A,cos_B,sin_A,sin_B,X,Y;
  double T0,TII,TT0,TT1;
  double R,S, PP,bang;
  double sin_TT1, sin_TT0, square=0.;
  double llai;



  llai = 0.;
  ang  = 5.*PI180;
  dang = 2.*ang;

  // Midpoints for azimuth intervals
  for (j=1; j<=16; j++)
    aden[j] = 0.0625;

  /* 3.1415/16 =0.1963 */
  for (j=1; j<=17; j++)
    {
      K         = 2*j - 3;
      TT[j]     = 0.1963*((double) K);
      sin_TT[j] = sin(TT[j]);
    }

  for (j=1; j<=16; j++)
    {
      del_TT[j]  = TT[j+1] - TT[j];
      del_sin[j] = sin_TT[j+1] - sin_TT[j];
    }

  for (KK = 1; KK <= 9; KK++)
    {
      bang = ang;
      // COMPUTE G FUNCTION FOR EACH LAYER IN THE CANOPY
      llai = 0.0;
      for (II=jtot; II>=1; II--)
	{
	  /*
	    'top of the canopy is jtot. Its cumulative LAI starts with 0.
	  */
	  llai += prof.dLAIdz[II];
	  // CALCULATE PROBABILITY FREQUENCY DISTRIBUTION, BDENS
	  FREQ(llai);
	  // LEMEUR DEFINES BDENS AS DELTA F/(PI/N), WHERE HERE N=9
	  PPP = 0.0;
	  for (I = 1; I <= 9; I++)
	    {
	      aang = ((I - 1) * 10.0 + 5.0) * PI180;
	      cos_B = cos(bang);
	      sin_B = sin(bang);
	      cos_A=cos(aang);
	      sin_A=sin(aang);
	      X = cos_A * sin_B;
	      Y = sin_A * cos_B;

	      if ((aang - bang) <= 0.0)
		{ /* if ab { */
		  for (IN = 1; IN <= 16; IN++)
		    {
		      /* for IN { */
		      PGF[IN] = X * del_TT[IN] + Y * del_sin[IN];
		    } /* for IN } */
		  goto LOOPOUT;
		} /* if ab } */
	      else
		{ /* else ab { */
		  T0 = 1.0 + X /Y;
		  TII = 1.0 - X / Y;

		  if (T0/TII > 0)
		    square = sqrt(T0/TII);
		  else
		    puts("bad sqrt in ggfun\n");

		  TT0 = 2.0 * atan(square);
		  sin_TT0 = sin(TT0);
		  TT1 = PI2 - TT0;
		  sin_TT1 = sin(TT1);

		  for (IN = 1,IN1=2; IN <= 16;IN1++, IN++)
		    { /* for IN { */
		      if ((TT[IN1] - TT0) <= 0)
			{ /* if 1a { */
			  PGF[IN] = X * del_TT[IN] + Y * del_sin[IN];
			  continue;
			} /* if 1a } */
		      else
			{ /* else 1a { */
			  if ((TT[IN1] - TT1) <= 0)
			    { /* if 1 { */
			      if ((TT0 - TT[IN]) <= 0)
				{ /* if 2 { */
				  PGF[IN] = -X * del_TT[IN] - Y * del_sin[IN];
				  continue;
				} /* if 2 } */
			      else
				{ /* else 2 { */
				  R = X * (TT0 - TT[IN]) + Y * (sin_TT0 - sin_TT[IN]);
				  S = X * (TT[IN1] - TT0) + Y * (sin_TT[IN1] - sin_TT0);
				  PGF[IN] = R - S;
				  continue;
				} /* else 2 } */
			    } /* if 1 } */
			  else
			    { /* else 1 { */
			      if ((TT1 - TT[IN]) <= 0.0)
				{ /* if 3 { */
				  PGF[IN] = X * del_TT[IN] + Y * del_sin[IN];
				  continue;
				} /* if 3 } */
			      else
				{ /* else 3 { */
				  R = X * (TT1 - TT[IN]) + Y * (sin_TT1 - sin_TT[IN]);
				  S = X * (TT[IN1] - TT1) + Y * (sin_TT[IN1] - sin_TT1);
				  PGF[IN] = S - R;
				} /* else 3 } */
			    } /* else 1 } */
			} /* else 1a } */
		    } /* else ab } */
		} /* next IN } */

	    LOOPOUT:

	      // Compute the integrated leaf orientation function, Gfun
	      PP = 0.0;
	      for (IN = 1; IN <= 16; IN++)
		PP += (PGF[IN] * aden[IN]);
	      PPP += (PP * canopy.bdens[I] * PI9);
	    } // next I
	  prof.Gfunc_sky[II][KK] = PPP;
	} // next IJ
      ang += dang;
    } // NEXT KK

  return;
}


// ----------------------------------------------------
void GFUNC()
{

  /*
    ------------------------------------------------------------
    This subroutine computes the G function according to the
    algorithms of Lemeur(1973). This progrom computes G for a given
    sun angle. G changes with height due to change leaf angles
    ------------------------------------------------------------
  */

  int IJ,IN, K, I,II, IN1;

  double aden[19], TT[19], pgg[19];
  double sin_TT[19], del_TT[19],del_sin[19];
  double PPP, PP, aang;
  double cos_A,cos_B,sin_A,sin_B,X,Y, sin_TT0, sin_TT1;
  double T0,TII,TT0,TT1;
  double R,S,square=0.;

  double llai = 0.;



  // Midpoint of azimuthal intervals

  for (IN=1; IN <= 17; IN++){
    K = 2 * IN - 3;
    TT[IN] = PI / 16.0 * K;
    sin_TT[IN] = sin(TT[IN]);
  }

  for (IN=1,IN1=2; IN <= 16; IN++,IN1++){
    del_TT[IN] = TT[IN1] - TT[IN];
    del_sin[IN]=sin_TT[IN1]-sin_TT[IN];
  }

  for (I=1; I <= 16; I++)
    aden[I] = .0625;


  // Compute the G function for each layer

  for (IJ = 1; IJ <= jtot; IJ++){
    II = jtot - IJ + 1;

    /* need LAI from LAIZ */

    llai += prof.dLAIdz[II];


    // Calculate the leaf angle probabilty distribution, bdens


    FREQ(llai);

    // Lemeur defines bdens as delta F/(PI/N), WHERE HERE N=9

    PPP = 0.0;

    for (I = 1; I <= 9; I++){
      aang = ((I - 1.0) * 10.0 + 5.0) * PI180;

      cos_A = cos(aang);
      cos_B = cos(solar.beta_rad);
      sin_A = sin(aang);
      sin_B = sin(solar.beta_rad);

      X = cos_A * sin_B;
      Y = sin_A * cos_B;

      if ((aang - solar.beta_rad) <= 0.0){
        for (IN = 1; IN <= 16; IN++){
          pgg[IN] = X * del_TT[IN] + Y * del_sin[IN];
        }
        goto OUTGF;
      }
      else{
        T0 = (1.0 + X/Y);
        TII = (1.0 - X/Y);

        if (T0/TII > 0)
          square=sqrt(T0/TII);
        else
          puts("bad T0/TII \n");


        TT0 = 2.0 * atan(square);
        TT1 = 2.0 * PI - TT0;
        sin_TT0=sin(TT0);
        sin_TT1=sin(TT1);

        for (IN = 1,IN1=2; IN <= 16; IN++,IN1++){
          if ((TT[IN1] - TT0) <= 0.0){
            pgg[IN] = X * del_TT[IN] + Y *del_sin[IN];
          }
          else{
            if ((TT[IN1] - TT1) <= 0.0){
              if ((TT0 - TT[IN]) <= 0.0){
                pgg[IN] = -X * del_TT[IN] - Y * del_sin[IN];
              }
              else{
                R = X * (TT0 - TT[IN]) + Y * (sin_TT0 - sin_TT[IN]);
                S = X * (TT[IN + 1] - TT0) + Y * (sin_TT[IN + 1] - sin_TT0);
                pgg[IN] = R - S;
              }
            }
            else{
              if ((TT1 - TT[IN]) <= 0.0){
                pgg[IN] = X * del_TT[IN] + Y * del_sin[IN];
              }
              else{
                R = X * (TT1 - TT[IN]) + Y * (sin_TT1 - sin_TT[IN]);
                S = X * (TT[IN1] - TT1) + Y * (sin_TT[IN1] - sin_TT1);
                pgg[IN] = S - R;
              }
            }
          }
        } // next IN
      }
    OUTGF:

      // Compute the integrated leaf orientation function

      PP = 0.0;
      for (IN = 1;IN <= 16; IN++)
        PP += (pgg[IN] * aden[IN]);


      PPP += (PP * canopy.bdens[I] * 9. / PI);

    } // next I

    prof.Gfunc_solar[II] = PPP;

    if (prof.Gfunc_solar[II] < 0.01)
      prof.Gfunc_solar[II] = 0.01;

    if (prof.Gfunc_solar[II] > 1.0)
      prof.Gfunc_solar[II] = 1.0;
  } // next IJ

  return;
}


// ----------------------------------------------------
void ENERGY_AND_CARBON_FLUXES()
{
  /*

  The ENERGY_AND_CARBON_FLUXES routine to computes coupled fluxes
  of energy, water and CO2 exchange, as well as leaf temperature. Computataions
  are performed for each layer in the canopy and on the sunlit and shaded fractions.

  Analytical solution for leaf energy balance and leaf temperature is used. The program
  is derived from work by Paw U (1986) and was corrected for errors with a re-derivation
  of the equations. The quadratic form of the solution is used, rather than the quartic
  version that Paw U prefers.

  Estimates of leaf temperature are used to adjust kinetic coefficients for enzyme kinetics,
  respiration and photosynthesis, as well as the saturation vapor pressure at the leaf surface.

  The Analytical solution of the coupled set of equations for photosynthesis and stomatal
  conductance by Baldocchi (1994, Tree Physiology) is used. This equation is a solution to
  a cubic form of the photosynthesis equation. The photosynthesis algorithms are from the
  model of Farquhar. Stomatal conductance is based on the algorithms of Ball-
  Berry and Collatz et al., which couple gs to A

  Layer 1 is the soil, Layer 40 is top of the canopy

  */

  int j;
  double Tair_K_filtered; // temporary absolute air temperature
  double T_sfc_K, T_sfc_C; // surface temperatures in Kelvin and Centigrade
  double H_sun, loutsun, Rn_sun, A_sun; // energy fluxes on sun leaves
  double H_shade, loutsh, Rn_shade, A_shade; // energy fluxes on shaded leaves
  double LE_leaf, LE_wet, H_leaf, lout_leaf;
  double wj_leaf=0., wc_leaf=0., surface_rh, surface_vpd;
  double rs_sun, rs_shade, A_mg, GPP, resp, internal_CO2, surface_CO2, chloroplast_CO2;
  double csca, cica, ccca;
  double fact_rs_sun, fact_rs_shd;

  for (j=1; j<=jtot; j++)
    {
      // zero summing values
      H_sun = 0.;
      LE_leaf = 0.;
      LE_wet = 0.;
      Rn_sun = 0.;
      loutsun = 0.;
      rs_sun = prof.sun_rs_filter[j];
      prof.sun_rs_save[j] = rs_sun;
      prof.sun_gs[j] = 1./rs_sun;
      A_sun = 0.;
      A_mg = 0.;
      GPP  = 0.;
      resp = 0.;
      internal_CO2 = prof.co2_air_filter[j];
      surface_CO2 = prof.co2_air_filter[j];
      chloroplast_CO2 = prof.co2_air_filter[j];
      csca = 1.;
      cica = 1.;
      ccca = 1.;
      surface_rh = prof.rhov_air_filter[j][1]*(prof.tair_filter[j]+TN0)*Rw
	/(ES(prof.tair_filter[j]+TN0)*100.);
      surface_vpd = (1.-surface_rh)*ES(prof.tair_filter[j]+TN0)*100.;
      /*
        First compute energy balance of sunlit leaf, then
        repeat and compute algorithm for shaded leaves.

        Remember layer is two-sided so include upper and lower streams
        are considered.

        KC is the convective transfer coeff. (W M-2 K-1). A factor
        of 2 is applied since convective heat transfer occurs
        on both sides of leaves.

        To calculate the total leaf resistance we must combine stomatal
        and the leaf boundary layer resistance. Many crops are amphistomatous
        so KE must be multiplied by 2. Deciduous forest, on the other hand
        has stomata on one side of the leaf.
      */

      Tair_K_filtered = prof.tair_filter[j] + TN0; // absolute air temperature

      // Initialize surface temperature with air temperature

      T_sfc_K = Tair_K_filtered;

      // Energy balance on sunlit leaves

      // update latent heat with new temperature during each call of this routine

      fact.latent = LAMBDA(Tair_K_filtered);
      /* if (j==1 || j==40) { printf("EC01 %i %20.17f %20.14f %20.14f\n", j, surface_rh, T_sfc_K, fact.latent); }; */
      if (solar.prob_beam[j] > 0.)
        {
          // Compute the resistances for heat and vapor transfer, rh and rv,
          // for each layer, s/m

          BOUNDARY_RESISTANCE(prof.ht[j], prof.sun_tleaf_filter[j], prof.cws[j][1], j);
	  /* if (j==1 || j==40) { printf("EC02 %20.17f %20.14f %20.14f\n", prof.ht[j], prof.sun_tleaf_filter[j], prof.cws[j][1]); }; */

          // compute energy balance of sunlit leaves
	  
          ENERGY_BALANCE(solar.rnet_sun[j], &T_sfc_K, Tair_K_filtered, prof.rhov_air_filter[j][1],
                         bound_layer_res.vapor, rs_sun, &LE_leaf, &LE_wet, &H_leaf, &lout_leaf,
                         prof.wet_coef_filter[j]);
	  /* if (j==1 || j==40) { printf("EC03 %20.17f %20.14f %20.14f\n", T_sfc_K, LE_leaf, LE_wet); }; */
	  /* if (j==1 || j==40) { printf("EC04 %20.14f %20.14f\n", H_leaf, lout_leaf); }; */

          // compute photosynthesis of sunlit leaves if leaves have emerged

          if (prof.dLAIdz[j] > pai/jtot)
            PHOTOSYNTHESIS(solar.quantum_sun[j], &rs_sun, prof.ht[j],
                           prof.co2_air_filter[j], T_sfc_K, LE_leaf, &A_mg, &GPP,
                           &resp, &internal_CO2, &surface_CO2, &chloroplast_CO2, &cica, &ccca,
                           &surface_rh, &surface_vpd,
                           &wj_leaf, &wc_leaf, j);

          // Assign values of function to the LE and H source/sink strengths

          T_sfc_C = T_sfc_K-TN0; // surface temperature, Centigrade
          H_sun = H_leaf; // sensible heat flux

	  prof.sun_LEstoma[j][1] = LE_leaf;    // latent heat flux from stomata
	  prof.sun_LEwet[j][1] = LE_wet;       // latent heat flux from wet leaf surface

          prof.sun_tleaf[j] = T_sfc_C;
          loutsun = lout_leaf; // long wave out
          Rn_sun = solar.rnet_sun[j] - lout_leaf; // net radiation

          A_sun = A_mg; // leaf photosynthesis, mg CO2 m-2 s-1

          prof.sun_resp[j] = resp; // respiration on sun leaves
          prof.sun_gs[j] = 1./rs_sun; // stomatal conductance
          prof.sun_rs[j] = rs_sun;

          fact_rs_sun = rugc * (prof.sun_tleaf[j]+TN0)/met.press_Pa; //conversion factor
          prof.sun_gs_mol[j] = 1./(prof.sun_rs[j]*fact_rs_sun*1.577); // convert to mol m-2 s-1 CO2

          prof.sun_wj[j] = wj_leaf;
          prof.sun_wc[j] = wc_leaf;

          prof.sun_GPP[j] = GPP; // micromolC m-2 s-1
          prof.sun_A[j] = A_sun*1000/mass_CO2; // micromolC m-2 s-1
          prof.sun_rbh[j] = bound_layer_res.heat;
          prof.sun_rbv[j] = bound_layer_res.vapor;
          prof.sun_rbco2[j] = bound_layer_res.co2;
          prof.sun_ci[j] = internal_CO2;
          prof.sun_cs[j] = surface_CO2;
          prof.sun_cc[j] = chloroplast_CO2;
          prof.sun_csca[j] = surface_CO2/prof.co2_air_filter[j];
          prof.sun_cica[j] = cica;
          prof.sun_ccca[j] = ccca;
          prof.sun_rh[j] = surface_rh; // relative humidity at leaf surface [0 to 1]
          prof.sun_vpd[j] = surface_vpd; // vapor pressure deficit at leaf surface [hPa]
        }

      // Energy balance on shaded leaves

      T_sfc_K = Tair_K_filtered;

      // initial value of stomatal resistance based on light

      rs_shade = prof.shd_rs_filter[j];
      prof.shd_rs_save[j] = rs_shade; // stomatal resistance, shaded leaf

      // boundary layer resistances on shaded leaves. With different
      // surface temperature, the convective effect may differ from that
      // computed on sunlit leaves
      BOUNDARY_RESISTANCE(prof.ht[j], prof.shd_tleaf_filter[j], prof.cws[j][1], j);
      /* if (j==1 || j==40) { printf("EC05 %20.17f %20.14f %20.14f\n", prof.ht[j], prof.shd_tleaf_filter[j], prof.cws[j][1]); }; */

      // Energy balance of shaded leaves
      ENERGY_BALANCE(solar.rnet_shd[j], &T_sfc_K, Tair_K_filtered, prof.rhov_air_filter[j][1],
                     bound_layer_res.vapor, rs_shade, &LE_leaf, &LE_wet, &H_leaf, &lout_leaf,
                     prof.wet_coef_filter[j]);
      /* if (j==1 || j==40) { printf("EC06.0 %20.17f %20.14f %20.14f\n", solar.rnet_shd[j], Tair_K_filtered, prof.rhov_air_filter[j][1]); }; */
      /* if (j==1 || j==40) { printf("EC06.1 %20.17f %20.14f %20.14f\n", bound_layer_res.vapor, rs_shade, prof.wet_coef_filter[j]); }; */
      /* if (j==1 || j==40) { printf("EC06 %20.17f %20.14f %20.14f\n", T_sfc_K, LE_leaf, LE_wet); }; */
      /* if (j==1 || j==40) { printf("EC07 %20.14f %20.14f\n", H_leaf, lout_leaf); }; */

      // compute photosynthesis and stomatal conductance of shaded leaves

      if (prof.dLAIdz[j] > pai/jtot)
        PHOTOSYNTHESIS(solar.quantum_shd[j], &rs_shade,prof.ht[j], prof.co2_air_filter[j],
                       T_sfc_K, LE_leaf, &A_mg, &GPP, &resp, &internal_CO2, &surface_CO2,
                       &chloroplast_CO2, &cica, &ccca,
                       &surface_rh, &surface_vpd, &wj_leaf,&wc_leaf, j);

      // re-assign variable names from functions output

      T_sfc_C = T_sfc_K-TN0; // surface temperature, C

      prof.shd_LEstoma[j][1] = LE_leaf;    // latent heat flux from stomata
      prof.shd_LEwet[j][1] = LE_wet;       // latent heat flux from wet leaf surface

      H_shade = H_leaf; // sensible heat flux density, shaded leaves, W m-2
      loutsh = lout_leaf; // long wave energy emissive flux density, shaded leaves, W m-2
      Rn_shade = solar.rnet_shd[j] - lout_leaf; // net radiation balance, shaded leaves, W m-2

      prof.shd_wj[j] = wj_leaf; // electron transport velocity, shaded leaves, micromol m-2 s-1
      prof.shd_wc[j] = wc_leaf; // carboxylation velocity, shaded leaves, micromol m-2 s-1

      A_shade = A_mg; // photosynthesis, shaded leaves, mgCO2 m-2 s-1

      prof.shd_rh[j] = surface_rh; // relative humidity at leaf surface [0 to 1]
      prof.shd_vpd[j] = surface_vpd; // vapor pressure deficit at leaf surface [hPa]

      // compute profiles

      prof.shd_resp[j] = resp;
      prof.shd_tleaf[j] = T_sfc_C;
      prof.shd_gs[j] = 1./rs_shade; // stomatal conductance, shaded leaf
      prof.shd_rs[j] = rs_shade; // stomatal resistance, shaded leaf

      fact_rs_shd = rugc * (prof.shd_tleaf[j]+TN0)/ met.press_Pa; //conversion factor
      prof.shd_gs_mol[j] = 1./(prof.shd_rs[j]*fact_rs_shd*1.577); // convert to mol m-2 s-1 CO2

      prof.shd_GPP[j] = GPP; // micromolC m-2 s-1
      prof.shd_A[j] = A_shade*1000/mass_CO2; // micromolC m-2 s-1
      prof.shd_rbh[j] = bound_layer_res.heat;
      prof.shd_rbv[j] = bound_layer_res.vapor;
      prof.shd_rbco2[j] = bound_layer_res.co2;
      prof.shd_ci[j]= internal_CO2;
      prof.shd_cs[j]= surface_CO2;
      prof.shd_cc[j] = chloroplast_CO2;
      prof.shd_csca[j] = surface_CO2/prof.co2_air_filter[j];
      prof.shd_cica[j] = cica;
      prof.shd_ccca[j] = ccca;

      // compute layer energy fluxes, weighted by leaf area and sun and shaded fractions

      prof.dLEdz_sun[j] = prof.dLAIdz[j] * solar.prob_beam[j]
	* (prof.sun_LEstoma[j][1]+prof.sun_LEwet[j][1]);
      prof.dLEdz_shd[j] = prof.dLAIdz[j] * solar.prob_shd[j]
	* (prof.shd_LEstoma[j][1]+prof.shd_LEwet[j][1]);
      prof.dLEdz[j][1] = prof.dLEdz_sun[j] + prof.dLEdz_shd[j];
      
      prof.dHdz[j] = prof.dLAIdz[j] * (solar.prob_beam[j] * H_sun
						  + solar.prob_shd[j] * H_shade);
      prof.dRNdz[j] = prof.dLAIdz[j] * (solar.prob_beam[j] * Rn_sun
						   + solar.prob_shd[j] * Rn_shade);
      prof.dLoutdz[j] = prof.dLAIdz[j] * (solar.prob_beam[j] * loutsun
						     + solar.prob_shd[j] * loutsh);

      // photosynthesis of the layer, prof.dPsdz has units mg m-3 s-1

      //prof.dPsdz[j] = prof.dLAIdz[j] * (A_sun * solar.prob_beam[j] + A_shade * solar.prob_shd[j]);

      // photosynthesis of layer, prof.dPsdz has units of micromoles m-2 s-1

      prof.dPsdz_sun[j] = prof.dLAIdz[j] * prof.sun_A[j] * solar.prob_beam[j];
      prof.dPsdz_shd[j] = prof.dLAIdz[j] * prof.shd_A[j] * solar.prob_shd[j];
      prof.dPsdz[j] = prof.dPsdz_sun[j] + prof.dPsdz_shd[j];

    // GPP of layer, prof.dGPPdz has units of micromoles m-2 s-1

      prof.dGPPdz_sun[j] = prof.dLAIdz[j] * prof.sun_GPP[j] * solar.prob_beam[j];
      prof.dGPPdz_shd[j] = prof.dLAIdz[j] * prof.shd_GPP[j] * solar.prob_shd[j];
      prof.dGPPdz[j] = prof.dGPPdz_sun[j] + prof.dGPPdz_shd[j];

      // respiration of the layer, micromol m-2 s-1

      prof.dRESPdz_sun[j] = prof.dLAIdz[j] * prof.sun_resp[j] * solar.prob_beam[j];
      prof.dRESPdz_shd[j] = prof.dLAIdz[j] * prof.shd_resp[j] * solar.prob_shd[j];
      prof.dRESPdz[j] = prof.dRESPdz_sun[j] + prof.dRESPdz_shd[j];

      // prof.dStomCondz has units of: m s-1
      // prof.dStomCondz has units of: mol m-2 s-1

      prof.dStomCondz_sun[j] = prof.dLAIdz[j] * solar.prob_beam[j]*prof.sun_gs[j];
      prof.dStomCondz_shd[j] = prof.dLAIdz[j] * solar.prob_shd[j]*prof.shd_gs[j];
      prof.dStomCondz[j] = prof.dStomCondz_sun[j] + prof.dStomCondz_shd[j];

      if ((prof.sun_LEwet[j][1]*solar.prob_beam[j]+
	   prof.shd_LEwet[j][1]*solar.prob_shd[j]) != 0.)
	prof.wet_coef[j]  = __max(__min(prof.cws[j][1] *
					       fact.latent/time_var.time_step/prof.dLAIdz[j]/
					       (prof.sun_LEwet[j][1]*solar.prob_beam[j]+
						prof.shd_LEwet[j][1]*solar.prob_shd[j])
					       , 1.), 0.);
      else
	prof.wet_coef[j] = 0.;
    } // next j

  return;
}


// ----------------------------------------------------
double SKY_IR(double T)
{
  // Infrared radiation from sky, W m-2, using algorithm from Norman
  double y;

  //y = sigma * pow(T,4) * ((1. - 0.261*exp(-.000777 * pow((TN0-T),2)))*solar.ratrad + 1.-solar.ratrad);
  y = sigma * T*T*T*T * ((1. - 0.261*exp(-.000777 * (TN0-T)*(TN0-T)))*solar.ratrad + 1.-solar.ratrad);

  return y;
}


// ----------------------------------------------------
void INIT_NEXT_STEP()
{
  int I, j, mc;

  // zero global arrays and variables

  // output
  output.c1                   = 0.;
  output.c2                   = 0.;
  output.c3                   = 0.;
  output.c4                   = 0.;
  output.c5                   = 0.;
  output.c6                   = 0.;
  output.c7                   = 0.;
  output.c8                   = 0.;
  output.c9                   = 0.;
  output.c10                  = 0.;
  output.netrad               = 0.;
  output.sumrn                = 0.;
  output.sumlout              = 0.;
  output.sumh                 = 0.;
  output.sumle                = 0.;
  output.can_ps_mol           = 0.;
  output.can_gpp              = 0.;
  output.c_transpiration_mole = 0.;
  output.canresp              = 0.;
  output.sumksi               = 0.;
  output.tleaf_mean           = 0.;
  output.tavg_sun             = 0.;
  output.tavg_shd             = 0.;
  output.ave_D13              = 0.;
  output.ave_D13_long         = 0.;
  output.ave_cica             = 0.;
  output.ave_gs               = 0.;
  output.ave_daC13            = 0.;
  output.diff_par             = 0.;
  output.rnet_soil            = 0.;
  for (I=1; I<=jtot; I++)
    output.dLAIdz[I]          = 0.;

  // flux
  flux.photosyn = 0.; // photosynthesis, um m-2 s-1
  flux.nee      = 0.; // net ecosystem exchange, umol m-2 s-1
  for (mc=1; mc<=wiso.nwater+1; mc++)
    {
      flux.c_evaporation[mc]        = 0.; // canopy evaporation from interception reservoir, mm m-2 s-1
      flux.c_transpiration[mc]      = 0.; // canopy transpiration through stomata, mm m-2 s-1
      flux.c_evapotranspiration[mc] = 0.; // canopy transpiration + evaporation, mm m-2 s-1
      flux.soilevap[mc]             = 0.; // soil evaporation, mm m-2 s-1
      flux.litterevap[mc]           = 0.; // litter evaporation, mm m-2 s-1
      flux.s_evap[mc]               = 0.; // soil+litter evaporation, mm m-2 s-1
      flux.evapotranspiration[mc]   = 0.; // canopy+soil+litter evaporationtranspiration, mm m-2 s-1
      flux.soilinfl[mc]             = 0.; // soil infiltration, mm m-2 s-1
      flux.surfrun[mc]              = 0.; // surface runoff = water flux lateral to surface, mm m-2 s-1
    }

  // fact
  fact.latent   = 0.; // latent heat of vaporization, J kg-1
  fact.heatcoef = 0.; // factor for sensible heat flux density
  fact.a_filt   = 0.; // filter coefficients
  fact.co2      = 0.; // CO2 factor, ma/mc * rhoa (mole m-3)

  // bole
  bole.factor           = 0.; // base respiration, micromoles m-2 s-1, data from Edwards
  bole.respiration_mole = 0.; // bole respiration, micromol m-2 s-1
  bole.respiration_mg   = 0.; // bole respiration, mg CO2 m-2 s-1
  bole.calc             = 0.; // calculation factor
  for (I=1; I<=jtot; I++)
    bole.layer[I]       = 0.; // bole pai per layer

  // solar
  for (I=1; I<=jtot1; I++)
    {
      solar.par_down[I]      = 0.;
      solar.par_up[I]        = 0.;
      solar.nir_dn[I]        = 0.;
      solar.nir_up[I]        = 0.;
      solar.ir_dn[I]         = 0.;
      solar.ir_up[I]         = 0.;
      solar.prob_beam[I]     = 0.;
      solar.prob_shd[I]      = 0.;
      solar.beam_flux_par[I] = 0.;
      solar.beam_flux_nir[I] = 0.;
      solar.rnet_sun[I]      = 0.;
      solar.rnet_shd[I]      = 0.;
      solar.nir_sun[I]       = 0.;
      solar.nir_shd[I]       = 0.;
      solar.par_shd[I]       = 0.;
      solar.par_sun[I]       = 0.;
      solar.quantum_sun[I]   = 0.;
      solar.quantum_shd[I]   = 0.;
    }

  // soil
  soil.T_Kelvin = 0.; // soil surface temperature in Kelvin
  soil.T_air = 0.; // air temperature above soil, C
  soil.Temp_ref = 0.; // reference soil temperature, annual mean, C
  soil.Temp_amp = 0.; // amplitude of soil temperature, C
  soil.rs = 0.; // soil resistance
  soil.rb = 0.; // soil boundary layer resistance
  soil.T_15cm = 0.; // soil temperature at 15 cm
  soil.lout = 0.; // longwave efflux from soil
  soil.evap = 0.; // total forest floor evaporation (soil+litter)
  soil.heat = 0.; // soil sensible heat flux density
  soil.rnet = 0.; // net radiation budget of the soil
  soil.gsoil = 0.; // soil heat flux density
  soil.respiration_mole = 0.; // soil respiration, micromol m-2 s-1
  soil.respiration_mg = 0.; // soil respiration, mg CO2 m-2 s-1
  soil.base_respiration = 0.; // base rate of soil respiration, micromol m-2 s-1
  soil.resp_13 = 0.; // respiration of 13C micromole m-2 s-1
  soil.SR_ref_temp = 0.; // refernce temperature for soil respiration, degC
  soil.rl = 0.; // litter resistance
  soil.litterevap_save = 0.; // save litter evaporation if it hits maximum
  soil.maxlitterevap = 0.; // maximum possible litter evaporation
  soil.c_litterevap = 0.; // factor regulating that litter evaporation < litter moisture
  soil.c_litterevap_save = 0.; // save factor regulating that litter evaporation < litter moisture
  soil.maxlitter = 0; // switch if litterevap was restricted to maximum litter evaporation 
  soil.latent = 0.; // latent heat coefficient before soil-air vapour pressure deficit
  soil.lecoef = 0.; // latent heat coefficient before soil-air vapour pressure deficit
  soil.rh_soil = 0.; // relative humidity in first soil layer air space
  soil.rh_litter = 0.; // relative humidity in litter layer air space
  soil.soil_mm = 0.; //water content in the soil, mm m-2
  soil.soil_mm_50 = 0.; //water content in the soil, mm m-2

  for (I=1; I<=soil.n_soil; I++)
    {
      soil.k_conductivity_soil[I] = 0.; // thermal conductivity of soil *dz
      soil.cp_soil[I] = 0.; // specific heat of soil, f(texture, moisture)
      soil.swp[I] = 0.; //soil water potential (kPa)
      soil.swp_mm[I] = 0.; //soil water potential (mm)
      soil.a[I] = 0.;
      soil.b[I] = 0.;
      soil.c[I] = 0.;
      soil.psi[I] = 0.; // soil water potential [Pa]
      soil.k_theta[I] = 0.; // moisture conductivity [mm s-1]
    }

  for (mc=1; mc<=wiso.nwater+1; mc++)
    {
      soil.qdrai[mc] = 0.; //sub-surface runoff (mm h2o s-1)
      soil.qinfl[mc] = 0.; //infiltration rate, mm h2o s-1
      soil.qtran[mc] = 0.; //plant transpiration, mm s-1
      soil.qseva[mc] = 0.; //plant transpiration, mm s-1
      soil.xylem[mc] = 0.; // xylem water isotopes
      soil.rxylem[mc] = 0.; // xylem water isotopic composition
    }

  for (mc=1; mc<=wiso.nwater+1; mc++)
    for (I=1; I<=soil.n_soil; I++)
      {
        soil.r[I][mc] = 0.;
        soil.dwat[I][mc] = 0.;
      }

  // prof
  for (I=1; I<=jtot; I++)
    {
      prof.sun_lai[I] = 0.;
      prof.shd_lai[I] = 0.;
      prof.dStomCondz[I] = 0.;
      prof.dStomCondz_sun[I] = 0.;
      prof.dStomCondz_shd[I] = 0.;
      prof.dPsdz[I] = 0.;
      prof.dPsdz_sun[I] = 0.;
      prof.dPsdz_shd[I] = 0.;
      prof.dGPPdz[I] = 0.;
      prof.dGPPdz_sun[I] = 0.;
      prof.dGPPdz_shd[I] = 0.;
      prof.dRNdz[I] = 0.;
      prof.dLoutdz[I] = 0.;
      prof.dRESPdz[I] = 0.;
      prof.dRESPdz_sun[I] = 0.;
      prof.dRESPdz_shd[I] = 0.;
      prof.dHdz[I] = 0.;
      bole.layer[I] = 0.;
      prof.sun_ci[I] = 0.;
      prof.shd_ci[I] = 0.;
      prof.sun_D13[I] = 0.;
      prof.shd_D13[I] = 0.;
      prof.shd_cica[I] = 0.;
      prof.sun_cica[I] = 0.;
      prof.source_co2[I] = 0.;
      prof.sun_wc[I] = 0.;
      prof.shd_wc[I] = 0.;
      prof.sun_wj[I] = 0.;
      prof.shd_wj[I] = 0.;
      prof.sour13co2[I] = 0.;
      prof.sun_D13_long[I] = 0.;
      prof.shd_D13_long[I] = 0.;
      prof.D13C_long[I] = 0.;
      prof.D13C_a[I] = 0.;
      prof.D13C_ab[I] = 0.;
      prof.D13C_asal[I] = 0.;
      prof.D13C_b[I] = 0.;
      prof.Gfunc_solar[I] = 0.; // leaf-sun direction cosine function
      prof.tleaf[I] = 0.;             // leaf temperature per layer (sun and shade)
      prof.isopreneflux[I] = 0.;      // isoprene flux per layer (sun and shade)
      prof.vcmax[I] = 0.;             // vcmax in per layer (temperature corrected)
      prof.jmax[I] = 0.;              // jmax per layer (temperature corrected)
      prof.rd[I] = 0.;                // rd per layer
      prof.vcmaxz[I] = 0.;            // vcmaxz per layer
      prof.wet_coef[I] = 0.;      // coefficient to regulate evaporation from wet leaf surface [0-1]
      prof.wet_coef_filter[I] = 0.;      // coefficient to regulate evaporation from wet leaf surface [0-1]
      prof.D13C[I] = 0.; // photosynthetic weighted discrimination, 13D
      prof.cs[I] = 0.; // CO2 concentration at leaf surface
      prof.ci[I] = 0.; // CO2 concentration inside stomatal cavity
      prof.cc[I] = 0.; // CO2 concentration at site of carboxlylation
      prof.d13Cplant[I] = 0.; // del 13C of the plant
      prof.Rplant_sun[I] = 0.; // ratio of discriminated 13C in sunlit leaves
      prof.Rplant_shd[I] = 0.; // ratio of discriminated 13C in shaded leaves
      prof.Rresp[I] = 0.; // discriminated 13C for later respiration, Ps weight
      prof.sun_frac[I] = 0.; // sun leaf fraction
      prof.sun_GPP[I] = 0.; // layer GPP flux for sun only (micromol mn-2 s-1)
      prof.sun_A[I] = 0.; // layer A flux for sun only (micromol mn-2 s-1)
      prof.sun_gs[I] = 0.; // stomatal conductance (m s-1)
      prof.sun_gs_mol[I] = 0.; // stomatal conductance of sun leaves mol m-2 s-1
      prof.sun_rbh[I] = 0.; // boundary layer resistance to heat (s/m)
      prof.sun_rbv[I] = 0.; // boundary layer resistance to H2O (s/m)
      prof.sun_rbco2[I] = 0.; // boundary layer resistance to CO2 (s/m)
      prof.sun_cs[I] = 0.; // Cs on sun leaves (CO2 mixing ratio on leaf surface
      prof.sun_cc[I] = 0.; // Cc (CO2 mixing ratio at site of carboxylation) on sun leaves
      prof.sun_D13_ab[I] = 0.; // del 13C boundary resistence on shaded leaves
      prof.sun_D13_a[I] = 0.; // del 13C stomatal resistance on shaded leaves
      prof.sun_D13_asal[I] = 0.;// del 13C mesophyll resistance on shaded leaves
      prof.sun_D13_b[I] = 0.; // del 13C rubisco on shaded leaves
      prof.sun_ccca[I] = 0.; // Cc/Ca on sunlit leaves
      prof.sun_resp[I] = 0.; // respiration
      prof.sun_isopreneflux[I] = 0.; // isoprene flux per layer for sunleaves
      prof.iso_sun[I] = 0.; // isoprene flux per leaf area in the sun
      prof.sun_vpd[I] = 0.; // vapor pressure deficit at leaf surface [hPa]
      prof.sun_rh[I] = 0.; // relative humidity at leaf surface [hPa]
      prof.shd_frac[I] = 0.; // shade leaf fraction
      prof.shd_GPP[I] = 0.; // photosynthesis (GPP) of shaded leaves
      prof.shd_A[I] = 0.; // photosynthesis of shaded leaves
      prof.shd_gs[I] = 0.; // stomatal conductance of shade leaves
      prof.shd_gs_mol[I] = 0.; // stomatal conductance of shade leaves mol m-2 s-1
      prof.shd_rbh[I] = 0.; // boundary layer resistance for heat on shade leaves
      prof.shd_rbv[I] = 0.; // boundary layer resistance for vapor on shade leaves
      prof.shd_rbco2[I] = 0.; // boundary layer resistance for CO2 on shade leaves
      prof.shd_cs[I] = 0.; // Cs on shade leaves (CO2 mixing ratio on leaf surface
      prof.shd_cc[I] = 0.; // Cc (CO2 mixing ratio at site of carboxylation) on sun leaves
      prof.shd_D13_ab[I] = 0.; // del 13C boundary resistence on shaded leaves
      prof.shd_D13_a[I] = 0.; // del 13C stomatal resistance on shaded leaves
      prof.shd_D13_asal[I] = 0.;// del 13C mesophyll resistance on shaded leaves
      prof.shd_D13_b[I] = 0.; // del 13C rubisco on shaded leaves
      prof.shd_ccca[I] = 0.; // Cc/Ca ratio on shaded leaves
      prof.shd_resp[I] = 0.; // respiration
      prof.shd_isopreneflux[I] = 0.; // isoprene flux per layer for shade leaves
      prof.iso_shd[I] = 0.; // isoprene flux per leaf area in the shade
      prof.shd_vpd[I] = 0.; // vapor pressure deficit at leaf surface [hPa]
      prof.shd_rh[I] = 0.; // relative humidity at leaf surface [hPa]
      prof.sun_wi[I] = 0.;
      prof.shd_wi[I] = 0.;
      prof.wa[I] = 0.;
      prof.sun_h[I] = 0.;
      prof.shd_h[I] = 0.;
      prof.rs_fact[I] = 0.;
      prof.sun_gross[I] = 0.;
      prof.shd_gross[I] = 0.;
      prof.sun_LEstoma_new[I] = 0.;
      prof.shd_LEstoma_new[I] = 0.;
      prof.sun_LEstoma_save[I] = 0.;
      prof.shd_LEstoma_save[I] = 0.;
    }

  for (mc=1; mc<=wiso.nwater+1; mc++)
    for (I=1; I<=jtot1; I++)
      prof.throughfall[I][mc] = 0.;     // throughfall in layer [mm], one extra layer for input as precipitation

  for (I=1; I<=jtot3; I++)
    {
      prof.d13C[I] = 0.;
      prof.c13cnc[I] = 0.;
      prof.c12cnc[I] = 0.;
      prof.u[I] = 0.; // wind speed (m/s)
      prof.vpd_air[I] = 0.; // vapor pressure deficit[kPa]
      prof.recycle[I] = 0.; // fraction of recycled CO2, after Yakir and Sternberg
    }

  for (mc=1; mc<=wiso.nwater+1; mc++)
    for (I=1; I<=jtot; I++)
      {
        prof.dLEdz[I][mc] = 0.;
	prof.dLEdz_sun[I] = 0.;
	prof.dLEdz_shd[I] = 0.;
        prof.sun_LEstoma[I][mc] = 0.;            // LE flux through stomata [W m-2 leaf]
        prof.sun_LEwet[I][mc] = 0.;         // LE flux from wet leaves [W m-2 leaf]
        prof.shd_LEstoma[I][mc] = 0.;            // LE flux through stomata [W m-2 leaf]
        prof.shd_LEwet[I][mc] = 0.;         // LE flux from wet leaves [W m-2 leaf]
	prof.sun_alpha_k[I][mc] = 0.;
	prof.shd_alpha_k[I][mc] = 0.;
	prof.sun_alpha_equ[I][mc] = 0.;
	prof.shd_alpha_equ[I][mc] = 0.;
	prof.sun_peclet[I][mc] = 0.;
	prof.shd_peclet[I][mc] = 0.;
	prof.sun_fem[I][mc] = 0.;
	prof.shd_fem[I][mc] = 0.;
	prof.sun_craig[I][mc] = 0.;
	prof.shd_craig[I][mc] = 0.;
	prof.sun_leafwater[I][mc] = 0.;
	prof.shd_leafwater[I][mc] = 0.;
	prof.sun_trans_rtrans[I][mc] = 0.;
	prof.shd_trans_rtrans[I][mc] = 0.;
	prof.sun_rtrans[I][mc] = 0.;
	prof.shd_rtrans[I][mc] = 0.;
      }

  // Copy former time step to current time step

  // Leaf water isotopes at evaporating site
  for (mc=1; mc<=wiso.nwater+1; mc++)
    for (j=1; j<=jtot; j++)
      {
	prof.sun_leafwater_e_old[j][mc] = prof.sun_leafwater_e[j][mc];
	prof.shd_leafwater_e_old[j][mc] = prof.shd_leafwater_e[j][mc];
	prof.sun_leafwater_e[j][mc] = 0.;
	prof.shd_leafwater_e[j][mc] = 0.;
      }
  
  // Do NOT zero the following variables!
  // They have to keep their values from the preceding time step
  /*
      // Surface temperature
      soil.tsfc
      soil.tsfc_filter
      // Interception reservoir
      for (mc=1; mc<=wiso.nwater+1; mc++)
	for (j=1; j<=jtot; j++)
	  prof.cws[j][mc]
      // z/L & sensible heat flux
      met.H
      met.H_filter
      // soil moisture
      for (mc=1; mc<=wiso.nwater+1; mc++)
	for (j=1; j<=soil.n_soil; j++)
          soil.theta[j][mc]
      soil.soil_mm_root
      // soil temperature
      for (j=0; j<=soil.n_soil+1; j++)
	{
	  soil.T_soil[j]
	  soil.T_soil_filter[j]
	}
      // litter moisture
      for (mc=1; mc<=wiso.nwater+1; mc++) // not executed if no water isotopes
        soil.theta_l[mc]
      // litter temperature
      soil.T_l
      soil.T_l_filter
      // profiles
      for (j=1; j<=jtot; j++)
	{
	  prof.sun_tleaf[j]
	  prof.sun_tleaf_filter[j]
	  prof.shd_tleaf[j]
	  prof.shd_tleaf_filter[j]
	}
      for (j=1; j<=jtot; j++)
	{
	  prof.sun_rs[j]
	  prof.sun_rs_filter[j]
	  prof.shd_rs[j]
	  prof.shd_rs_filter[j]
	}
      for (j=1; j<=jtot3; j++)
	{
	  prof.tair[j]
	  prof.tair_filter[j]
	  prof.co2_air[j]
	  prof.co2_air_filter[j]
	}
      for (mc=1; mc<=wiso.nwater+1; mc++)
	for (j=1; j<=jtot3; j++)
	  {
	    prof.rhov_air[j][mc]
	    prof.rhov_air_filter[j][mc]
	  }
      if (set_switch.d13c == 1)
	{
	  for (j=1; j<=jtot3; j++)
	    {
	      prof.R13_12_air[j]
	      prof.d13Cair[j]
	    }
	}
      for (mc=2; mc<=wiso.nwater+1; mc++)
	for (j=1; j<=jtot3; j++)
	  prof.rvapour[j][mc]
      // evapotranspiration
      soil.soilevap
      soil.soilevap_filter
      soil.litterevap
      soil.litterevap_filter
  */
  return;
}


// ----------------------------------------------------
void ENERGY_BALANCE(double qrad, double *tsfckpt, double taa, double rhovva, double rvsfc,
                    double stomsfc, double *lept, double *lewet, double *H_leafpt, double *lout_leafpt,
                    double wet_coef)
{
  /*
    ENERGY BALANCE COMPUTATION


    A revised version of the quadratic solution to the leaf energy balance relationship is used.

    Paw U, KT. 1987. J. Thermal Biology. 3: 227-233


    H is sensible heat flux density on the basis of both sides of a leaf
    J m-2 s-1 (W m-2). Note KC includes a factor of 2 here for heat flux
    because it occurs from both sides of a leaf.

    LE is latent heat flux density through stomata. For hypostomatous leaves (n_stomata_sides = 1) LE
    occurs on only one side. For amphistomatous leaves (n_stomata_sides = 2) LE occurs
    on both sides.
    LEw is latent heat flux density from wet leave surfaces. It occurs on both sides.
  */

  double est, ea, tkta;
  double tk2, tk3, tk4;
  double dest, d2est;
  double lecoef, hcoef, hcoef2, product;
  //double bcoef, ccoef, repeat,acoef, acoeff, le2;
  double atlf, btlf, ctlf,vpd_leaf,llout;
  double ke;

  tkta = taa;               //[K]
  est = ES(tkta) * 100.;   // converts es(T) from mb to Pa

  // ea = RHOA * TAA * 1000 / 2.165

  ea = rhovva * taa * Rw; // vapor pressure above leaf

  // Vapor pressure deficit, Pa

  vpd_leaf = est - ea;

  // Slope of the vapor pressure-temperature curve, Pa/C
  // evaluate as function of Tk

  dest = DESDT(tkta);
  /* printf("EB01 %20.14f %20.14f %20.14f\n", tkta, fact.latent, dest); */

  // Second derivative of the vapor pressure-temperature curve, Pa/C
  // Evaluate as function of Tk

  d2est = DES2DT(tkta);
  /* printf("EB02 %20.14f\n", d2est); */

  // Compute products of air temperature, K

  tk2 = tkta * tkta;
  tk3 = tk2 * tkta;
  tk4 = tk3 * tkta;

  // Longwave emission at air temperature, W m-2

  llout = epsigma2 * tk4;
  /* printf("EB03 %20.14f %20.14f %20.14f\n", epsigma2, tk4, llout); */

  /*
    Check if leaves have stomata on one or both sides
    Cuticle resistance is included in STOM.

    hypostomatous  n_stomata_sides = 1
    amphistomatous  n_stomata_sides = 2

    read in from parameter file in set_parameter()
  */
  /*
    Interception led to wet leaf surface (both sides of leaves), which can contribute
    to evaporation from leaf surfaces. Therefore we have to include an evaporaton term
    without stomata resistance (only boundary layer resistance) We introduce
    wet_coef[0-1] to make sure that LE from wet surface does not exceed the amount
    of water on the leaves.
  */

  //  Coefficient for latent heat flux

  ke = 1./ (rvsfc + stomsfc) + wet_coef/rvsfc ;

  lecoef = met.air_density * 0.622 * fact.latent * ke / met.press_Pa;

  // Coefficients for sensible heat flux

  hcoef = met.air_density*cp/bound_layer_res.heat;
  hcoef2 = 2. * hcoef;

  /* now LE is not directly calculated with the quadratic solution anymore,
     but we first solve for leaf temperature with a quadratic equation. Then we use
     the outcome of leaf temperature to calculate LE

     // The quadratic coefficients for the a LE^2 + b LE +c =0

     repeat = hcoef + epsigma4 * tk3;

     acoeff = lecoef * d2est / (2. * repeat);
     acoef = acoeff / 4.;

     bcoef = -(repeat) - lecoef * dest / 2. + acoeff * (-qrad / 2. + llout);
     ccoef = repeat * lecoef * vpd_leaf + lecoef * dest * (qrad / 2. - llout) + acoeff * ((qrad * qrad) / 4. + llout * llout - qrad * llout);

     // LE1 = (-BCOEF + (BCOEF ^ 2 - 4 * ACOEF * CCOEF) ^ .5) / (2 * ACOEF)

     product = bcoef * bcoef - 4. * acoef * ccoef;

     // LE2 = (-BCOEF - (BCOEF * BCOEF - 4 * acoef * CCOEF) ^ .5) / (2. * acoef)

     le2= (-bcoef - sqrt(product)) / (2. * acoef);

     *lept=le2; // need to pass pointer out of subroutine
  */

  // solve for Ts using quadratic solution

  // coefficients to the quadratic solution

  atlf = epsigma12 * tk2 + n_stomata_sides * d2est * lecoef / 2.;

  btlf = epsigma8 * tk3 + hcoef2 + n_stomata_sides * lecoef * dest;

  ctlf = -qrad + llout + n_stomata_sides * lecoef * vpd_leaf;

  product = btlf * btlf - 4 * atlf * ctlf;
  /* printf("EB04 %20.14f %20.14f\n", lecoef, hcoef); */

  if (product >= 0.)
    *tsfckpt = tkta + (-btlf + sqrt(product)) / (2. * atlf); // [K]
  else
    *tsfckpt = tkta; // [K]


  if (*tsfckpt < 230. || *tsfckpt > 335.)
    *tsfckpt = tkta; // [K]

  // long wave emission of energy
  //*lout_leafpt = epsigma2 * pow(*tsfckpt, 4);
  *lout_leafpt = llout
    + epsigma8 * tkta*tkta*tkta * (*tsfckpt-tkta)
    + epsigma12 * tkta*tkta * (*tsfckpt-tkta)*(*tsfckpt-tkta);

  // H is sensible heat flux
  *H_leafpt = hcoef2 * (*tsfckpt-tkta);

  // lept is latent heat flux through stomata
  // ToDo for isotopes // transpiration from second Taylor expansion
  //*lept = n_stomata_sides * met.air_density * 0.622 * fact.latent
  //  / (met.press_Pa * (rvsfc + stomsfc)) * (vpd_leaf + dest*(*tsfckpt-tkta) + d2est/2.*(*tsfckpt-tkta)*(*tsfckpt-tkta));
  *lept = n_stomata_sides * met.air_density * 0.622 * fact.latent
    / (met.press_Pa * (rvsfc + stomsfc)) * (ES(*tsfckpt)*100.-ea);

  if (set_switch.no_neg_water_flux == 1)
    *lept = __max(*lept, 0.);

  // lewet is latent heat flux from wet leaves
  // ToDo for isotopes // soil evaporation from second Taylor expansion
  //*lewet = wet_coef * n_stomata_sides*met.air_density * 0.622 * fact.latent
  //  / (met.press_Pa * rvsfc) * (vpd_leaf + dest*(*tsfckpt-tkta) + d2est/2.*(*tsfckpt-tkta)*(*tsfckpt-tkta));
  *lewet = wet_coef * n_stomata_sides*met.air_density * 0.622 * fact.latent
    / (met.press_Pa * rvsfc) * (ES(*tsfckpt)*100.-ea);
  *lewet = __max(*lewet, 0.);

  return;
}


// ----------------------------------------------------
double LAMBDA(double tak)
{ // Latent heat of Vaporisation, J kg-1
  double y;

  y = 3149000. - 2370. * tak;
  // add heat of fusion for melting ice
  if (tak < TN0)
    y +=333;

  return y;
}



// ----------------------------------------------------
double ES(double t)
{
  // saturation vapor pressure function (mb)
  // T is temperature in Kelvin
  double y=0.;

  if (t > 0.)
    {
      y = 54.8781919 - 6790.4985 / t - 5.02808 * log(t);
      y = exp(y);
    }
  else
    puts("bad es calc");

  return y;
}


// ----------------------------------------------------
double DESDT (double t)
{
  // first derivative of es with respect to tk
  // routine needs to convert es(t) (mb) to Pa
  double y;

  y = ES(t)*100. * fact.latent*18. / (rgc1000 *t*t);

  return y;
}


// ----------------------------------------------------
double DES2DT (double t)
{
  // The second derivative of the saturation vapor pressure
  // temperature curve, using the polynomial equation of Paw U
  // a3en=1.675;
  // a4en=0.01408;
  // a5en=0.0005818;
  double y,tcel;

  tcel = t-TN0;

  y = 2.*a3en + 6.*a4en*tcel + 12.*a5en*tcel*tcel;

  return y;
}


// ----------------------------------------------------
void PHOTOSYNTHESIS(double Iphoton,double *rstompt, double zzz, double cca, double tlk,
                    double leleaf, double *A_mgpt, double *GPPpt, double *resppt, double *cipnt, double *cspnt,
                    double *ccpnt, double *cicapnt, double *cccapnt,
                    double *rh_leafpnt, double *vpd_leafpnt, double *wjpnt, double *wcpnt, int JJ)
{
  /*
    This program solves a cubic equation to calculate
    leaf photosynthesis. This cubic expression is derived from solving
    five simultaneous equations for A, PG, cs, CI and GS.
    Stomatal conductance is computed with the Ball-Berry model.
    The cubic derivation assumes that b', the intercept of the Ball-Berry
    stomatal conductance model, is non-zero.

    Gs = k A rh/cs + b'

    We also found that the solution for A can be obtained by a quadratic equation
    when Gs is constant or b' is zero.

    The derivation is published in:

    Baldocchi, D.D. 1994. An analytical solution for coupled leaf photosynthesis
    and stomatal conductance models. Tree Physiology 14: 1069-1079.

    -----------------------------------------------------------------------

    A Biochemical Model of C3 Photosynthesis

    After Farquhar, von Caemmerer and Berry (1980) Planta.
    149: 78-90.

    The original program was modified to incorporate functions and parameters
    derived from gas exchange experiments of Harley, who paramertized Vc and J in
    terms of optimal temperature, rather than some reference temperature, eg 25C.

    Program calculates leaf photosynthesis from biochemical parameters

    rd25 - Dark respiration at 25 degrees C (umol m-2 s-1)
    tlk - leaf temperature, Kelvin
    jmax - optimal rate of electron transport
    vcopt - maximum rate of RuBP Carboxylase/oxygenase
    iphoton - incident photosynthetically active photon flux (mmols m-2 s-1)

    note: Harley parameterized the model on the basis of incident PAR

    gs - stomatal conductance (mols m-2 s-1), typically 0.01-0.20
    pstat-station pressure, bars
    aphoto - net photosynthesis (umol m-2 s-1)
    ps - gross photosynthesis (umol m-2 s-1)
    aps - net photosynthesis (mg m-2 s-1)
    aphoto (umol m-2 s-1)

    --------------------------------------------------

    iphoton is radiation incident on leaves

    The temperature dependency of the kinetic properties of
    RUBISCO are compensated for using the Arrhenius and
    Boltzmann equations. From biochemistry, one observes that
    at moderate temperatures enzyme kinetic rates increase
    with temperature. At extreme temperatures enzyme
    denaturization occurs and rates must decrease.

    Arrhenius Eq.

    f(T)=f(tk_25) exp(tk -298)eact/(298 R tk)), where eact is the
    activation energy.

    Boltzmann distribution

    F(T)=tboltz()

    Define terms for calculation of gross photosynthesis, PG

    PG is a function of the minimum of RuBP saturated rate of
    carboxylation, Wc, and the RuBP limited rate of carboxylation, Wj.
    Wj is limiting when light is low and electron transport, which
    re-generates RuBP, is limiting. Wc is limiting when plenty of RuBP is
    available compared to the CO2 that is needed for carboxylation.

    Both equations take the form:

    PG-photorespiration= (a CI-a d)/(e CI + b)

    PG-photorespiration=min[Wj,Wc] (1-gamma/Ci)

    Wc=Vcmax Ci/(Ci + Kc(1+O2/Ko))

    Wj=J Ci/(4 Ci + 8 gamma)

    Ps kinetic coefficients from Harley at WBW.

    Gamma is the CO2 compensation point

    Jan 14, 1999 Updated the cubic solutions for photosynthesis. There are
    times when the restriction that R^2 < Q^3 is violated. I therefore need
    alternative algorithms to solve for the correct root.

    ===============================================================
  */

  double tprime25, bc, ttemp, gammac;
  double jmax, vcmax, jmaxz=0., vcmaxz=0., cs, ci, cc, vc25z;
  double kct, ko, tau;
  double rd, rdz;
  double rb_mole, gb_mole, dd, b8_dd;
  double rh_leaf, vpd_leaf, es_leaf, k_rh, ci_guess;
  double j_photon, alpha_ps=0., bbeta=0., ggamma=0.;
  double denom, Pcube=0., Qcube=0., Rcube=0.;
  double P2, P3, Q, R;
  double root1, root2;
  double root3, arg_U, ang_L;
  double aphoto=0., gpp=0., j_sucrose, wj=0.;
  double gs_leaf_mole=0., gs_co2, gs_m_s;
  double ps_1,delta_1, Aquad1, Bquad1, Cquad1;
  double theta_ps=0., wc=0., b_ps, a_ps, e_ps, psguess;
  double sqrprod=0., product; // delday;
  double rt;
  double gm;
  double beta_ps=0., gamma_ps=0., delta_ps=0., epsilon_ps=0., zeta_ps=0., eta_ps=0.;
  double theta_soil;
  double g0_local=0., a1_local=0., D0_local=0.;

  // double a_cubic, b_cubic, rootprod;
  double rr, qqq, minroot=0., maxroot=0., midroot=0.;
  double bprime_local=0., bprime16_local;

  rt = rugc * tlk; // product of universal gas constant and abs temperature

  tprime25 = tlk - tk_25; // temperature difference

  ttemp = exp((skin * tlk - hkin) / rt) + 1.0; // denominator term

  if (Iphoton < 1.) Iphoton = 0.;

  // KC and KO are solely a function of the Arrhenius Eq.

  kct = TEMP_FUNC(kc25, ekc, tprime25, tk_25, tlk); // åµbar
  ko = TEMP_FUNC(ko25, eko, tprime25, tk_25, tlk); // mbar
  tau = TEMP_FUNC(tau25, ektau, tprime25, tk_25, tlk); // dimensonless

  bc = kct * (1.0 + o2 / ko); // åµbar*(1+mbar/mbar) = åµbar

  /*
    gammac is the CO2 compensation point due to photorespiration, umol mol-1
    Recalculate gammac with the new temperature dependent KO and KC
    coefficients

    gammac = .5 * O2*1000/TAU
  */

  gammac = 500.0 * o2 / tau; //0.5*mmol mol-1 * 100/tau = åµmol mol-1 = ppm

  /*
    temperature corrections for Jmax and Vcmax

    Scale jmopt and VCOPT with a surrogate for leaf nitrogen
    specific leaf weight (Gutschick and Weigel).

    normalized leaf wt is 1 at top of canopy and is 0.35
    at forest floor. Leaf weight scales linearly with height
    and so does jmopt and vcmax
    zoverh=0.65/HT=zh65

  */

  if (extra_nate == 1)
    {
      // for Nate McDowell's juniper site, no scaling of Vcmax with height and LAI
      jmaxz = jmopt[JJ];
      vcmaxz = vcopt[JJ];
    }
  else
    {
      // time before leaf out
      if (time_var.days < time_var.leafout)
	{
	  jmaxz = 0.;
	  vcmaxz = 0.;
	}
      
      // spring, increase Ps capacity with leaf expansion as a function of leaf area changes
      if (time_var.days >= time_var.leafout && time_var.days < time_var.leaffull)
	{
	  jmaxz = jmopt[JJ] * (zh65 * zzz + .35) * time_var.lai/lai;
	  vcmaxz = vcopt[JJ] * (zh65 * zzz + .35) * time_var.lai/lai;
	}
      
      // growing season, full Ps capacity (note newer data by Wilson et al shows more
      // dynamics      
      if (time_var.days >= time_var.leaffull && time_var.days < time_var.leaffall)
	{
	  jmaxz = jmopt[JJ] * (zh65 * zzz + .35);
	  vcmaxz = vcopt[JJ] * (zh65 * zzz + .35);
	}
      
      // gradual decline in fall
      if (time_var.days >= time_var.leaffall && time_var.days <= time_var.leaffallcomplete)
	{
	  //delday=1-(time_var.days-270)/30;
	  jmaxz = jmopt[JJ] * (zh65 * zzz + .35) * time_var.lai/lai;
	  vcmaxz = vcopt[JJ] * (zh65 * zzz + .35) * time_var.lai/lai;
	}
      
      if (time_var.days > time_var.leaffallcomplete)
	{
	  jmaxz = 0.;
	  vcmaxz = 0.;
	}
    }

  prof.vcmaxz[JJ] = vcmaxz;

  /*
    Scale rd with height via vcmax and apply temperature
    correction for dark respiration
  */
  // get vcmax at 25 deg C and then take fraction of it as dark respiration (rd_vc)
  vc25z = TBOLTZ(vcmaxz, evc, toptvc, TN0+25.);
  rdz = vc25z * rd_vc[JJ];

  //? from Harley 1995, sun leaves: Rd(25 deg C)/Vcmax(Topt)=0.34/73=0.004657?
  // but if we use Rd(25 deg C)/Vcmax(25 deg C)=0.34/34=0.01 (Harley 1995 data)
  // collatz 1991 gives rd=0.015*vcmax
  // Farqhuar 1980 gives rd=0.011*vcmax


  // reduce respiration by 50% in light according to Amthor (might be less, see Pinelli and Loreto, 2003)
  if (Iphoton > 10) rdz *= 0.5; //changed to 40% reduction

  // apply temperature correction for rd at 25 deg C to leaf level temperature
  rd = TEMP_FUNC(rdz, erd, tprime25, tk_25, tlk);

  prof.rd[JJ] = rd; //store rd in gobal structure

  // Apply temperature correction to JMAX and vcmax

  jmax = TBOLTZ(jmaxz, ejm, toptjm, tlk);
  vcmax = TBOLTZ(vcmaxz, evc, toptvc, tlk);

  prof.jmax[JJ] = jmax; //store jmax in gobal structure
  prof.vcmax[JJ] = vcmax; //store vcmax in gobal structure

  /*
    Compute the leaf boundary layer resistance

    gb_mole leaf boundary layer conductance for CO2 exchange,
    mol m-2 s-1

    RB has units of s/m, convert to mol-1 m2 s1 to be
    consistant with R.

    rb_mole = RBCO2 * .0224 * 1.01 * tlk / (met.pstat * TN0)
  */

  rb_mole = bound_layer_res.co2 * tlk * (met.pstat273); // met.pstat273 = rugc / (100000 * met.press_bars)

  gb_mole = 1. / rb_mole;

  dd = gammac;
  b8_dd = 8. * dd;


  /*
    **************************************

    aphoto = PG - rd, net photosynthesis is the difference
    between gross photosynthesis and dark respiration. Note
    photorespiration is already factored into PG.

    ****************************************

    coefficients for Ball-Berry stomatal conductance model
    
    Gs = k A rh/cs + b'

    rh is relative humidity, which comes from a coupled
    leaf energy balance model
  */

  rh_leaf = SFC_VPD(tlk, zzz, leleaf); // calculate relative humidity at leaf surface

  es_leaf = ES(tlk); // saturation vapor pressure at leaf temperature
  vpd_leaf = es_leaf - rh_leaf*es_leaf; // calculate vapor pressure deficit at leaf surface


  if (soil.camillo == 0)
    {
      // sets response of stomata conductance to soil water stress
      // calculate relative plant available water
      // normalise total soil water by total soil water at field capacity (-33 kPa)
      theta_soil = __min(__max((soil.soil_mm_root-soil.soil_mm_1500_root)
			       /(soil.soil_mm_33_root-soil.soil_mm_1500_root), 0.), 1.);

      output.c9 = theta_soil;

      if (theta_soil > sfc_res.fthreshold)
        sfc_res.fdrought = 1.;
      else // prevent instability by setting minimum to 0
        sfc_res.fdrought = __min(__max(1.-sfc_res.fslope*(sfc_res.fthreshold-theta_soil), 0.), 1.);
    }
  else
    sfc_res.fdrought = 1.;


  // set Ball-Berry stomatal factor
  k_rh = sfc_res.fdrought * rh_leaf * kball[JJ]; // combine product of rh and K ball-berry times drought response factor

  /*
    Gs from Ball-Berry is for water vapor. It must be divided
    by the ratio of the molecular diffusivities to be valid
    for A
  */
  k_rh = k_rh / 1.577; // adjust the coefficient for the diffusion of CO2 rather than H2O


  // parameter for Leuning mesophyll model
  if (extra_nate == 1)
    { // !!! ASK Alex about 1.577 etc. !!!
      g0_local = g0[JJ]; //empirical coefficient, intercept, converted from water to CO2
      //      a1_local = sfc_res.fdrought*a1[JJ]; //slope of stomata function [-] converted from water to CO2 times drought response factor
      a1_local = a1[JJ]; //slope of stomata function [-] converted from water to CO2 times drought response factor
      D0_local = D0[JJ]; //empirical coefficient, intercept, converted from water to CO2
    }
  else
    {
      g0_local = g0[JJ]/1.577; //empirical coefficient, intercept, converted from water to CO2
      a1_local = sfc_res.fdrought*a1[JJ]/1.577; //slope of stomata function [-] converted from water to CO2 times drought response factor
      D0_local = D0[JJ]; //empirical coefficient, intercept, converted from water to CO2
    }

  gm = gm_vc*vcmax; //mesophyll conductance [mol m-2 s-1]
  prof.gm[JJ] = gm;

  ci_guess = cca * 0.7; // initial guess of internal CO2 to estimate Wc and Wj

  if (Iphoton < 1.) goto quad; // shortcut to dark

  /*
    Test for the minimum of Wc and Wj. Both have the form:

    W = (a ci - ad)/(e ci + b)

    after the minimum is chosen set a, b, e and d for the cubic solution.

    estimate of J according to Farquhar and von Cammerer (1981)


    J photon from Harley
  */

  if (jmax > 0)
    //j_photon = qalpha * Iphoton / sqrt(1. +(qalpha2 * Iphoton * Iphoton / (jmax * jmax))); // Michaelis Menten type (originally in CANOAK)
    // non rectangular hyperbola with curvature parameter (see Von Cammerer 2000)
    //    j_photon = (qalpha * Iphoton + jmax - sqrt(pow(qalpha*Iphoton+jmax,2)-4.*curvature*qalpha*jmax*Iphoton)) / (2.*curvature);
    j_photon = (qalpha * Iphoton + jmax - sqrt((qalpha*Iphoton+jmax)*(qalpha*Iphoton+jmax)-4.*curvature*qalpha*jmax*Iphoton)) / (2.*curvature);
  else
    j_photon = 0.;

  wj = j_photon * (ci_guess - dd) / (4. * ci_guess + b8_dd);

  wc = vcmax * (ci_guess - dd) / (ci_guess + bc);

  // frost and end of leaf photosynthesis and respiration

  if (time_var.days > time_var.leaffallcomplete) // old: 300
    {
      wj = 0.;
      j_photon = 0.;
      wc = 0.;
      rd = 0.;
    }

  if (wj < wc)
    {
      // for Harley and Farquhar type model for Wj
      psguess = wj;

      b_ps = b8_dd;
      a_ps = j_photon;
      e_ps = 4.;
    }
  else
    {
      psguess=wc;

      b_ps = bc;
      a_ps = vcmax;
      e_ps = 1.;
    }



  // cubic coefficients that are only dependent on CO2 levels

  // for the Ball Berry Farquhar model based on Baldocchis analytical solution
  if (set_switch.ball == 0)
    {
      bprime_local = bprime[JJ];
      bprime16_local = bprime[JJ]/1.577;
      alpha_ps = 1.0 + (bprime16_local / gb_mole) - k_rh;
      bbeta = cca * (gb_mole * k_rh - 2.0 * bprime16_local - gb_mole);
      ggamma = cca * cca * gb_mole * bprime16_local;
      theta_ps = gb_mole * k_rh - bprime16_local;
    }
  else  // for the Leuning Farquhar model based on Knohls analytical solution
    {
      alpha_ps = 1/gm + 1/gb_mole;
      beta_ps = (D0_local+vpd_leaf)*gb_mole*(cca-gammac);
      gamma_ps = a1_local*gb_mole*vpd_leaf;
      delta_ps = b_ps+e_ps*cca;
      epsilon_ps = a_ps-e_ps*rd;
      zeta_ps = D0_local+vpd_leaf;
      eta_ps = a_ps*(gammac-cca);
      /*alpha_ps = (gb_mole*gm+g0*(gm+gb_mole))*(D0+vpd_leaf);
        beta_ps = gb_mole*gm*g0*(D0+vpd_leaf);
        gamma_ps = gb_mole*D0*a1;
        delta_ps = a_ps*dd-a_ps*cca+rd*e_ps*cca+rd*e_ps;
        epsilon_ps = gb_mole*gm*gamma_ps-beta_ps;
        zeta_ps = gm*gamma_ps+gb_mole*gamma_ps-alpha_ps;
      */
      bprime_local = bprime[JJ];
      bprime16_local = bprime[JJ]/1.577;
    }

  /*
    if wj or wc are less than rd then A would probably be less than zero. This would yield a
    negative stomatal conductance. In this case, assume gs equals the cuticular value. This
    assumptions yields a quadratic rather than cubic solution for A
  */

  if (wj <= rd) goto quad;

  if (wc <= rd) goto quad;

  /*
    cubic solution: A^3 + p A^2 + q A + r = 0
  */

  // for the Ball Berry Farquhar model based on Baldocchis analytical solution
  if (set_switch.ball == 0)
    {
      denom = e_ps * alpha_ps;

      Pcube = (e_ps * bbeta + b_ps * theta_ps - a_ps * alpha_ps + e_ps * rd * alpha_ps);
      Pcube = Pcube/denom;

      Qcube = (e_ps * ggamma + (b_ps * ggamma / cca) - a_ps * bbeta + a_ps * dd * theta_ps + e_ps * rd * bbeta + rd * b_ps * theta_ps);
      Qcube = Qcube/denom;

      Rcube = (-a_ps * ggamma + a_ps * dd * (ggamma / cca) + e_ps * rd * ggamma + rd * b_ps * ggamma / cca);
      Rcube = Rcube/denom;
    }
  else  // for the Leuning Farquhar model based on Knohls analytical solution
    {

      denom = (-alpha_ps*e_ps*gamma_ps+alpha_ps*e_ps*g0_local*zeta_ps+e_ps*zeta_ps);

      Pcube = alpha_ps*epsilon_ps*gamma_ps-alpha_ps*epsilon_ps*g0_local*zeta_ps-alpha_ps*e_ps*g0_local*beta_ps-zeta_ps*epsilon_ps-beta_ps*e_ps+delta_ps*gamma_ps-delta_ps*zeta_ps*g0_local;
      Pcube = Pcube/denom;

      Qcube = alpha_ps*epsilon_ps*g0_local*beta_ps+beta_ps*epsilon_ps+delta_ps*beta_ps*g0_local+eta_ps*gamma_ps-eta_ps*zeta_ps*g0_local+delta_ps*gamma_ps*rd-delta_ps*zeta_ps*g0_local*rd;
      Qcube = Qcube/denom;

      Rcube = eta_ps*beta_ps*g0_local+delta_ps*beta_ps*g0_local*rd;
      Rcube = Rcube/denom;
      /*
        denom = e_ps * (-zeta_ps);

        Pcube = (a_ps-e_ps*rd)*zeta_ps-e_ps*alpha_ps*cca*gb_mole+(b_ps+e_ps*cca)*epsilon_ps;
        Pcube = Pcube/denom;

        Qcube = (b_ps+e_ps*cca)*beta_ps*cca*gb_mole+delta_ps*epsilon_ps;
        Qcube = Qcube/denom;

        Rcube = delta_ps*beta_ps*cca*gb_mole;
        Rcube = Rcube/denom;
      */
    }


  // Use solution from Numerical Recipes from Press

  P2 = Pcube * Pcube;
  P3 = P2 * Pcube;
  Q = (P2 - 3.0 * Qcube) / 9.0;
  R = (2.0 * P3 - 9.0 * Pcube * Qcube + 27.0 * Rcube) / 54.0;

  /*
    Test = Q ^ 3 - R ^ 2
    if test >= O then all roots are real
  */

  rr = R*R;
  qqq = Q*Q*Q;

  // real roots

  arg_U = R / sqrt(qqq);

  ang_L = acos(arg_U);

  root1 = -2.0 * sqrt(Q) * cos(ang_L / 3.0) - Pcube / 3.0;
  root2 = -2.0 * sqrt(Q) * cos((ang_L + PI2) / 3.0) - Pcube / 3.0;
  root3 = -2.0 * sqrt(Q) * cos((ang_L - PI2) / 3.0) - Pcube / 3.0;

  // rank roots #1,#2 and #3 according to the minimum, intermediate and maximum
  // value

  if (root1 <= root2 && root1 <= root3)
    {
      minroot=root1;
      if (root2 <= root3)
        {
          midroot=root2;
          maxroot=root3;
        }
      else
        {
          midroot=root3;
          maxroot=root2;
        }
    }


  if (root2 <= root1 && root2 <= root3)
    {
      minroot=root2;
      if (root1 <= root3)
        {
          midroot=root1;
          maxroot=root3;
        }
      else
        {
          midroot=root3;
          maxroot=root1;
        }
    }


  if (root3 <= root1 && root3 <= root2)
    {
      minroot=root3;
      if (root1 < root2)
        {
          midroot=root1;
          maxroot=root2;
        }
      else
        {
          midroot=root2;
          maxroot=root1;
        }

    } // end of the loop for real roots

  // find out where roots plop down relative to the x-y axis

  if (minroot > 0. && midroot > 0. && maxroot > 0.) aphoto = minroot;

  if (minroot < 0. && midroot < 0. && maxroot > 0.) aphoto = maxroot;

  if (minroot < 0. && midroot > 0. && maxroot > 0.) aphoto = midroot;

  /*
    Here A = x - p / 3, allowing the cubic expression to be expressed
    as: x^3 + ax + b = 0
  */

  // aphoto=root3; // back to original assumption

  /*
    also test for sucrose limitation of photosynthesis, as suggested by
    Collatz. Js=Vmax/2
  */
  j_sucrose = vcmax / 2. - rd;

  if (j_sucrose < aphoto) aphoto = j_sucrose;

  cs = cca - aphoto / gb_mole;

  if (cs > 3.*cca) cs = input.co2air;

  gpp = aphoto + rd;

  /* Stomatal conductance for water vapor

  forest are hypostomatous.
  Hence we don't divide the total resistance
  by 2 since transfer is going on only one side of a leaf.
  */

  if (set_switch.ball == 0)
    gs_leaf_mole = (k_rh*1.577 * aphoto / cs) + bprime_local;
  else
    gs_leaf_mole = a1_local*1.577*aphoto/((cs-dd)*(1+vpd_leaf/D0_local)) + g0_local*1.577;


  // convert Gs from vapor to CO2 diffusion coefficient based on diffusivities in Massman (1998)
  gs_co2 = gs_leaf_mole / 1.577;

  if (aphoto < 0.) puts("aphoto<0 should not be here");

  ci = cs - aphoto / gs_co2;
  cc = ci - aphoto / gm;

  if (set_switch.ball == 0)
    {
      wj = j_photon * (ci - dd) / (4. * ci + b8_dd);
      wc = vcmax * (ci - dd) / (ci + bc);
    }
  else
    {
      wj = j_photon * (cc - dd) / (4. * cc + b8_dd);
      wc = vcmax * (cc - dd) / (cc + bc);
    }

  /* stomatal conductance is mol m-2 s-1
     convert back to resistance (s/m) for energy balance routine
  */

  gs_m_s = gs_leaf_mole * tlk * met.pstat273;

  // need point to pass rstom out of subroutine

  *rstompt = 1.0 / gs_m_s;


  //   // to compute ci, Gs must be in terms for CO2 transfer
  //   ci = cs - aphoto / gs_co2;

  //   // compute cc with mesophyll conductance
  //   cc = ci - aphoto / gm;

  /*
    if A < 0 then gs should go to cuticular value and recalculate A
    using quadratic solution
  */

  if (aphoto <= 0.) goto quad;

  goto OUTDAT;

  // if aphoto < 0 set stomatal conductance to cuticle value
 quad:

  gs_leaf_mole = bprime_local;
  gs_co2 = gs_leaf_mole / 1.577;

  /*
    stomatal conductance is mol m-2 s-1
    convert back to resistance (s/m) for energy balance routine
  */
  gs_m_s = gs_leaf_mole * tlk * (met.pstat273);

  // need point to pass rstom out of subroutine as a pointer

  *rstompt = 1.0 / gs_m_s;

  /*
    a quadratic solution of A is derived if gs=ax, but a cubic form occurs
    if gs =ax + b. Use quadratic case when A is less than zero because gs will be
    negative, which is nonsense

  */
  if (Iphoton < 1.) // shortcut dark
    {
      gpp = 0.;
      aphoto = -rd;
    }
  else
    {
      ps_1 = cca * gb_mole * gs_co2;
      delta_1 = gs_co2 + gb_mole;
      denom = gb_mole * gs_co2;

      Aquad1 = delta_1 * e_ps;
      Bquad1 = -ps_1 * e_ps - a_ps * delta_1 + e_ps * rd * delta_1 - b_ps * denom;
      Cquad1 = a_ps * ps_1 - a_ps * dd * denom - e_ps * rd * ps_1 - rd * b_ps * denom;

      product = Bquad1 * Bquad1 - 4.0 * Aquad1 * Cquad1;

      if (product >= 0.)
	{
	  sqrprod= sqrt(product);
	  aphoto = (-Bquad1 - sqrprod) / (2.0 * Aquad1);
	  /*
	    Tests suggest that APHOTO2 is the correct photosynthetic root when
	    light is zero because root 2, not root 1 yields the dark respiration
	    value rd.
	  */
	}
      else
	{
	  aphoto = 0.;
	}
      // correct for gpp>0 - should only be numerical adjustments
      gpp = __max(aphoto+rd, 0.);
      aphoto = gpp-rd;
    }

  cs = cca - aphoto / gb_mole;
  ci = cs - aphoto / gs_co2;
  cc = ci - aphoto / gm;


 OUTDAT:

  if (gpp < 0.) puts("gpp<0 should not be here");
  /*
    compute photosynthesis with units of mg m-2 s-1 and pass out as pointers
    A_mg = APHOTO * 44 / 1000
  */
  *A_mgpt = aphoto * mass_CO2/1000.;
  *GPPpt = gpp;
  *resppt = rd;
  *cipnt = ci;
  *cspnt = cs;
  *ccpnt = cc;
  *cicapnt = ci/cca;
  *cccapnt = cc/cca;
  *wcpnt = wc;
  *wjpnt = wj;
  *rh_leafpnt = rh_leaf;
  *vpd_leafpnt = vpd_leaf;

  return;
}


// ----------------------------------------------------
double GAMMAF(double x)
{
  // gamma function

  double y=0., gam=0.;

  //  gam= (1.0 / (12.0 * x)) + (1.0 / (288.0 * x*x)) - (139.0 / (51840.0 * pow(x,3.0)));
  gam = (1.0 / (12.0 * x)) + (1.0 / (288.0 * x*x)) - (139.0 / (51840.0 * x*x*x));
  gam += 1.;

  if (x > 0.)
    y = sqrt(2.0 * PI / x) * pow(x,x) * exp(-x) * gam;
  else
    puts("gamma\n");

  return y;
}


// ----------------------------------------------------
double TBOLTZ (double rate, double eakin, double topt, double tl)
{
  // Boltzmann temperature distribution for photosynthesis

  double y, dtlopt,prodt,numm,denom;

  dtlopt = tl - topt;
  prodt = rugc * topt * tl;
  numm = rate * hkin * exp(eakin * (dtlopt) / (prodt));
  denom = hkin - eakin * (1.0 - exp(hkin * (dtlopt) / (prodt)));
  y = numm / denom;
  return y;
}


// ----------------------------------------------------
double INV_TBOLTZ (double rate, double eakin, double topt, double tl)
{

  // calculates the inverse of Boltzmann temperature distribution for photosynthesis
  // used to get Vcmax and Jmax at Toptimum from Vcmax and Jmax at 25 deg C

  double y, dtlopt,prodt,numm,denom;

  dtlopt = tl - topt;
  prodt = rugc * topt * tl;
  denom = hkin * exp(eakin * (dtlopt) / (prodt));
  numm = hkin - eakin * (1.0 - exp(hkin * (dtlopt) / (prodt)));
  y = rate * (numm / denom);
  return y;
}


// ----------------------------------------------------
double TEMP_FUNC(double rate,double eact,double tprime,double tref, double t_lk)
{
  // Arhennius temperature function
  double y;
  y = rate * exp(tprime * eact / (tref * rugc*t_lk));
  return y;
}


// ----------------------------------------------------
void ANGLE()
{

  // ANGLE computes solar elevation angles,

  // This subroutine is based on algorithms in Walraven. 1978. Solar Energy. 20: 393-397

  double theta_angle,G,EL,EPS,sin_el,A1,A2,RA;
  double delyr,leap_yr,T_local,time_1980,leaf_yr_4, delyr4;
  double day_savings_time, day_local;
  double S,HS,phi_lat_radians,value,declination_ang=0.,ST,SSAS;
  double E_ang=0., zenith, elev_ang_deg, cos_zenith;

  double radd = 0.0174532925;
  double twopi = 6.2831853072;
    
  delyr = time_var.year - 1980.0;
  delyr4=delyr/4.0;
  leap_yr=fmod(delyr4,4.0);
  day_savings_time=0.0;

  // Daylight Savings Time, Dasvtm =1
  // Standard time, Dasvtm= 0

  T_local = time_var.local_time;

  day_local = time_var.days;
  time_1980 = delyr * 365.0 + leap_yr + (day_local - 1) + T_local / 24.0;

  leaf_yr_4=leap_yr*4.0;

  if (delyr == leaf_yr_4) time_1980 -= 1.0;

  if (delyr < 0.0)
    if (delyr < leaf_yr_4 || delyr > leaf_yr_4) time_1980 -= 1.0;

  theta_angle = (360.0 * time_1980 / 365.25) * radd;
  G = -.031272 - 4.53963E-7 * time_1980 + theta_angle;
  EL = 4.900968 + 3.6747E-7 * time_1980 + (.033434 - 2.3E-9 * time_1980) * sin(G) + .000349 * sin(2. * G) + theta_angle;
  EPS = .40914 - 6.2149E-9 * time_1980;
  sin_el = sin(EL);
  A1 = sin_el * cos(EPS);
  A2 = cos(EL);

  RA = atan(A1/A2);

  /* for ATAN2

  RA = atan2(A1,A2);
  if (RA < 0)
  RA=RA+twopi;

  */

  /*
    The original program was in FORTRAN. It used the ATAN2 function.

    In C we must find the correct quadrant and evaluate the arctan
    correctly.

    Note ATAN2 is -PI TO PI, while ATN is from PI/2 TO -PI/2

  */

  // QUAD II, TAN theta_angle < 0

  if (A1 > 0)
    if (A2 <= 0.)
      RA += PI;

  // QUAD III, TAN theta_angle > 0 /

  if (A1 <= 0.)
    if (A2 <= 0.)
      RA += PI;

  value = sin_el * sin(EPS);

  if (1.-value*value >= 0.)
    declination_ang = atan(value/ sqrt(1. - value * value));
  else
    puts(" bad declination_ang\n");

  // declination_ang=asin(value)

  ST = 1.759335 + twopi * (time_1980 / 365.25 - delyr) + 3.694E-7 * time_1980;

  if (ST >= twopi) ST = ST - twopi;

  S = ST - longitude * radd + 1.0027379 * (zone - day_savings_time + T_local) * 15. * radd;

  if (S >= twopi) S = S - twopi;

  HS = RA - S;
  phi_lat_radians = latitude * radd;

  // DIRECTION COSINE

  SSAS = (sin(phi_lat_radians) * sin(declination_ang) + cos(phi_lat_radians) * cos(declination_ang) * cos(HS));

  if (1.-SSAS*SSAS >= 0.)
    E_ang = atan(sqrt(1. - (SSAS * SSAS))/ SSAS);
  else
    puts(" bad SSAS \n");

  if (SSAS < 0) E_ang=E_ang+PI;

  // E=asin(SSAS);

  if (E_ang < 0) E_ang=PI/2.;

  zenith = E_ang / radd;
  elev_ang_deg = 90. - zenith;
  solar.beta_rad = elev_ang_deg * radd;
  solar.sine_beta = sin(solar.beta_rad);
  cos_zenith = cos(zenith);
  solar.beta_deg = solar.beta_rad / PI180;

  return;
}


// ----------------------------------------------------
void CONC(double *source, double *cncc, double cref, double soilflux, double factor)
{


  // Subroutine to compute scalar concentrations from source
  // estimates and the Lagrangian dispersion matrix


  double sumcc[sze3], cc[sze3];
  double disper, ustfact, disperzl, soilbnd;

  int i, j;


  // Compute concentration profiles from Dispersion matrix


  ustfact = ustar_ref/ met.ustar; // factor to adjust Dij with alternative u* values
  // Note that disperion matrix was computed using u* = 0.405

  for ( i=1; i<=jtot3; i++)
    {
      sumcc[i] = 0.0;


      for (j=1; j<=jtot; j++)
        {

          /*
            CC is the differential concentration (Ci-Cref)


            Ci-Cref = SUM (Dij S DELZ), units mg m-3,mole m-3, or heat m-3

            S = dfluxdz/DELZ

            note delz values cancel

            scale dispersion matrix according to friction velocity
          */

          disper = ustfact * met.dispersion[i][j]; // units s/m

          // scale dispersion matrix according to Z/L


          // if (met.zl < 0)
          // disperzl = disper / (1.- .812 * met.zl);
          // else
          // disperzl=disper;

          // updated Dispersion matrix (Dec, 2002). New Tl and 1,000,000 particles Hainich

          if (met.zl < 0.)
            disperzl = disper * (0.2469* met.zl -0.4701)/(met.zl -0.3716); //changed
          else
            disperzl = disper;

          sumcc[i] += delz * disperzl * source[j];


        } // next j


      // scale dispersion matrix according to Z/L

      // case for j=1, soil level

      disper = ustfact * met.dispersion[i][1];


      if (met.zl < 0.)
        disperzl = disper * (0.2469* met.zl -0.4701)/(met.zl -0.3716); //changed
      else
        disperzl = disper;


      // add soil flux to the lowest boundary condition
      // convert to the units we need

      soilbnd = soilflux*disperzl/factor;

      cc[i] = sumcc[i]/factor+soilbnd;


    } // next i


  // Compute scalar profile below reference

  for (i=1; i<=jtot3; i++)
    cncc[i] = cc[i] + cref - cc[izref];

  return;
}


// ----------------------------------------------------
void IRFLUX()
{
  int j;
  double reflc_lay_IR;
  double Tk_sun_filt,Tk_shade_filt, IR_source_sun,IR_source_shade,IR_source;
  double SDN[sze], SUP[sze], ir_dn[sze], ir_up[sze];
  double emiss_IR_soil;
  double r2, r22;
  int i_check;

  /*
    This subroutine is adapted from:
    Norman, J.M. 1979. Modeling the complete crop canopy.
    Modification of the Aerial Environment of Crops.
    B. Barfield and J. Gerber, Eds. American Society of Agricultural Engineers, 249-280.

    Compute probability of penetration for diffuse radiation for each layer in the canopy.
    IR radiation is isotropic.
  */

  // IR down flux at top of canopy
  solar.ir_dn[jtot1] = SKY_IR(met.T_Kelvin);

  /*
    Integrated probability of diffuse sky radiation penetration
    EXPDIF[JJ] is computed in RAD

    compute IR radiative source flux as a function of
    leaf temperature weighted according to sunlit and shaded fractions

    source=ep*sigma*(laisun*tksun^4 + laish*tksh^4)
    remember energy balance is done on layers not levels.
    so level jtot+1 must use tl from layer jtot
  */
  for (j=jtot; j>=1; j--)
    {
      Tk_sun_filt     = prof.sun_tleaf_filter[j]+TN0;
      Tk_shade_filt   = prof.shd_tleaf_filter[j]+TN0;
      IR_source_sun   = solar.prob_beam[j] * Tk_sun_filt*Tk_sun_filt*Tk_sun_filt*Tk_sun_filt;
      IR_source_shade = solar.prob_shd[j]  * Tk_shade_filt*Tk_shade_filt*Tk_shade_filt*Tk_shade_filt;
      IR_source       = epsigma * (IR_source_sun + IR_source_shade);
      // Intercepted IR that is radiated up
      SUP[j+1]        = IR_source * (1.-solar.exxpdir[j]);
      // Intercepted IR that is radiated downward
      SDN[j]          = IR_source * (1.-solar.exxpdir[j]);
    }

  // First guess of up and down values

  // Downward IR radiation, sum of that from upper layer that is transmitted
  //   and the downward source generated in the upper layer.
  for (j=jtot; j>=1; j--)
    solar.ir_dn[j] = solar.exxpdir[j]   * solar.ir_dn[j+1] + SDN[j];
  // same for upward
  for (j=2; j<=jtot1; j++)
    solar.ir_up[j] = solar.exxpdir[j-1] * solar.ir_up[j-1] + SUP[j];
  // ground emission
  emiss_IR_soil         = epsoil*sigma
    * (soil.tsfc_filter+TN0)*(soil.tsfc_filter+TN0)
    * (soil.tsfc_filter+TN0)*(soil.tsfc_filter+TN0);
  SUP[1]                = solar.ir_dn[1] * (1.-epsoil);
  solar.ir_up[1] = emiss_IR_soil + SUP[1];

  i_check = 0;
  do
    { // save result
      for (j=1; j<=jtot1; j++)
	{
	  ir_up[j] = solar.ir_up[j];
	  ir_dn[j] = solar.ir_dn[j];
	}
      // go down
      for (j=jtot; j>=1; j--)
	{
          reflc_lay_IR          = (1.-solar.exxpdir[j]) * (1.-ep);
          solar.ir_dn[j] = solar.exxpdir[j] * solar.ir_dn[j+1]
	    + solar.ir_up[j] * reflc_lay_IR + SDN[j];
	}
      // go up
      SUP[1]                = solar.ir_dn[1] * (1.-epsoil);
      solar.ir_up[1] = emiss_IR_soil + SUP[1];
      for (j=2; j<=jtot1; j++)
        {
          reflc_lay_IR          = (1.-solar.exxpdir[j-1]) * (1.-ep);
          solar.ir_up[j] = reflc_lay_IR * solar.ir_dn[j]
	    + solar.ir_up[j-1] * solar.exxpdir[j-1] + SUP[j];
        }
      // check if stabilised
      r2  = 1. - pow(corr(jtot1, &ir_up[1], &solar.ir_up[1]), 2);
      r22 = 1. - pow(corr(jtot1, &ir_dn[1], &solar.ir_dn[1]), 2);
      if (i_check++ > 10)
	puts("bad irflux");
    }
  while (r2 > 1e-12 || r22 > 1e-12);

  return;
}


// ----------------------------------------------------
void DIFFUSE_DIRECT_RADIATION()
{
  double fAND, fir,fv;
  double rdir, rdvis,rsvis,wa;
  double ru, rsdir, rvt,rit, nirx;
  double xvalue,fvsb,fvd,fansb;
  /*
    This subroutine uses the Weiss-Norman ( 1985, Agric. forest Meteorol. 34: 205-213)
    routine to compute direct and diffuse PAR from total par

    fractions of NIR and PAR (visible)
  */
  fir = 0.54;
  fv = 0.46;

  ru = 98.5 / (101.3 * solar.sine_beta);
  /*
    visible direct PAR
  */

  rdvis = 600. * exp(-0.185 * ru) * solar.sine_beta;

  /*
    potential diffuse PAR
  */
  // rsvis = .4 * (600.0 - rdvis) * solar.sine_beta;
  rsvis = 0.4 * (600. * solar.sine_beta - rdvis); //new equation after Al Weiss, 1/14/05


  /*
    solar constant: 1320 W m-2

    water absorption in NIR for 10 mm precip water
  */

  wa = 1320. * 0.077 * pow((2.*ru), 0.3);

  /*

  direct beam NIR
  */
  rdir = (720. * exp(-0.06 * ru) - wa) * solar.sine_beta;

  if (rdir < 0.)
    rdir = 0.;

  /*
    potential diffuse NIR
  */

  //rsdir = .6 * (720 - wa - rdir) * solar.sine_beta;
  rsdir = .6 * (720 - wa - rdir/solar.sine_beta) * solar.sine_beta; //new equation after Al Weiss, 1/14/05

  if (rsdir < 0.)
    rsdir = 0.;

  rvt = rdvis + rsvis;
  rit = rdir + rsdir;

  if (rit < 0.1)
    rit = 0.1;

  if (rvt < 0.1)
    rvt = 0.1;

  solar.ratrad = input.rglobal / (rvt + rit);
  output.c10 = (rvt + rit);

  if (time_var.local_time >= 12. && time_var.local_time < 13.)
    solar.ratradnoon = solar.ratrad;

  /*
    ratio is the ratio between observed and potential radiation

    NIR flux density as a function of PAR

    since NIR is used in energy balance calculations
    convert it to W m-2: divide PAR by 4.6
  */


  nirx = input.rglobal - (input.parin / 4.6);

  /*
    ratio = (PARIN / 4.6 + NIRX) / (rvt + rit)
  */
  if (solar.ratrad > 0.89)
    solar.ratrad = 0.89;

  if (solar.ratrad < 0.22)
    solar.ratrad = 0.22;

  /*
    fraction PAR direct and diffuse
  */
  xvalue=(0.9-solar.ratrad)/0.70;

  fvsb = rdvis / rvt * (1. - pow(xvalue, 0.67));

  if (fvsb < 0.)
    fvsb = 0.;

  if (fvsb > 1.)
    fvsb = 1.;

  fvd = 1. - fvsb;

  /*
    note PAR has been entered in units of umol m-2 s-1
  */
  // use modeled direct and diffuse relationship
  solar.par_beam = fvsb * input.parin;
  solar.par_diffuse = fvd * input.parin;

  // use calculated PAR for Nate McDowell's Juniper data
  if (extra_nate == 0)
    { // use input data
      solar.par_beam = input.parin-input.pardif;
      solar.par_diffuse = input.pardif;
    }

  if (solar.par_beam <= 0.)
    {
      solar.par_beam = 0.;
      solar.par_diffuse = input.parin;
    }

  if (input.parin <= 0.001)
    {
      solar.par_beam = 0.001;
      solar.par_diffuse = 0.001;
    }

  xvalue=(0.9-solar.ratrad)/0.68;
  fansb = rdir / rit * (1. - pow(xvalue,0.67));

  if (fansb < 0.) fansb = 0.;

  if (fansb > 1.) fansb = 1.;

  fAND = 1. - fansb;

  /*
    NIR beam and diffuse flux densities
  */
  solar.nir_beam = fansb * nirx;
  solar.nir_diffuse = fAND * nirx;

  if (solar.nir_beam <= 0.)
    {
      solar.nir_beam = 0.;
      solar.nir_diffuse = nirx;
    }

  if (nirx <= 0.1)
    {
      solar.nir_beam = 0.1;
      solar.nir_diffuse = 0.1;
    }

  return;
}


// ----------------------------------------------------
void FREQ(double lflai)
{
  int I;

  double STD, MEAN, CONS;
  double THETA2,nuu,SUM,MU,FL1,MU1,nu1;
  double ANG,FL2,FL3;
  double leaf_angle_top, leaf_angle_bottom, leaf_angle_SD_top, leaf_angle_SD_bottom; 


  /*
    THIS PROGRAM USES THE BETA DISTRIBUTION
    TO COMPUTE THE PROBABILITY FREQUENCY
    DISTRIBUTION FOR A KNOWN MEAN LEAF INCLINATION ANGLE
    STARTING FROM THE TOP OF THE CANOPY, WHERE llai=0

    AFTER GOEL AND STREBEL (1984)

    Note: there is a typo in equation [8] in Goel and Strebel (1984)
    correct form is:

    nuu = (1. - <theta>^2 / (90. * <theta>)) / (<theta^2> / <theta>^2 - 1.);

    theta : leaf angle
    <theta> : mean of leaf angle

    Canopy type MEAN THETA2 STD
    Planophile 26.8 1058.6 18.5
    Erectophile 63.2 4341.4 18.5
    Plagiophile 45.0 2289.7 16.3
    Extremophile 45.0 3110.4 32.9
    Uniform 45.0 2700.0 26.0
    Sperical 57.3 3747.6 21.5


    // old
    if (lflai < 2.4)    //till here plagiophile
    MEAN = 45.;
    else
    MEAN = 78.99 - 14.84 * lflai + .244 * lflai * lflai;

    if(MEAN < 26)        //from here planophile
    MEAN=26;
    
    // MEAN=70.;            // mean leaf angle for eucalyptus from King 1997
    STD = 18;
  */

  leaf_angle_top = 55;        // Mean leaf inclination angle at canopy top 
  leaf_angle_bottom = 25;    // Mean leaf inclination angle at canopy bottom, old=25
  leaf_angle_SD_top = 21.5;        // Standard deviation of leaf inclination angle at canopy top 
  leaf_angle_SD_bottom = 18;    // Standard deviation of leaf inclination angle at canopy bottom, old =18

  MEAN = leaf_angle_top;
  STD = leaf_angle_SD_top;

  if (lflai < 2.4)
    {
      MEAN = leaf_angle_top;
      STD = leaf_angle_SD_top;
    }
  if ((lflai >= 2.4) && (lflai < 3.5))
    {
      MEAN = leaf_angle_top - (lflai-2.4)/(3.5-2.4)*(leaf_angle_top-leaf_angle_bottom);
      STD = leaf_angle_SD_top - (lflai-2.4)/(3.5-2.4)*(leaf_angle_SD_top-leaf_angle_SD_bottom);
    }
  if (lflai >= 3.5)
    {
      MEAN = leaf_angle_bottom;
      STD = leaf_angle_SD_bottom;
    }

  if (extra_nate == 1)
    {
      MEAN = 57.3;
      STD = 21.5;
      //MC !!! v2
      MEAN = 45.;
      STD = 18.;
      //MC !!! smallest netrad error
      MEAN = 35.;
      STD = 18.;
    }

  THETA2 = STD * STD + MEAN * MEAN; // = <theta^2> since Variance = STD*STD = <theta^2> - <theta>^2
  nuu = (1. - THETA2 / (90. * MEAN)) / (THETA2 / (MEAN * MEAN) - 1.);
  MU = nuu * ((90. / MEAN) - 1.);
  SUM = nuu + MU;

  FL1 = GAMMAF(SUM) / (GAMMAF(nuu) * GAMMAF(MU));
  MU1 = MU - 1.;
  nu1 = nuu - 1.;

  CONS = 1. / 9.;

  /*
    COMPUTE PROBABILITY DISTRIBUTION FOR 9 ANGLE CLASSES
    BETWEEN 5 AND 85 DEGREES, WITH INCREMENTS OF 10 DEGREES
  */
  for (I=1; I <= 9;I++)
    {
      ANG = (10. * I - 5.);
      FL2 = pow((1. - ANG / 90.), MU1);
      FL3 = pow((ANG / 90.), nu1);
      canopy.bdens[I] = CONS * FL1 * FL2 * FL3;
    }
  return;
}


// ----------------------------------------------------
double SFC_VPD (double tlk, double Z, double leleafpt)
{

  // this function computes the relative humidity at the leaf surface for
  // application in the Ball Berry Equation

  // latent heat flux, LE, is passed through the function, mol m-2 s-1
  // and it solves for the humidity at leaf surface

  int J;
  double y, rhov_sfc,e_sfc,vpd_sfc,rhum_leaf;
  double es_leaf;

  es_leaf = ES(tlk); // saturation vapor pressure at leaf temperature

  J = (int)(Z / delz); // layer number

  rhov_sfc = (leleafpt / (fact.latent)) * bound_layer_res.vapor + prof.rhov_air[J][1]; /* kg m-3 */

  e_sfc = rhov_sfc * tlk / 0.2165; // mb
  vpd_sfc = es_leaf - e_sfc; // mb
  rhum_leaf = 1. - vpd_sfc / es_leaf; // 0 to 1
  y = rhum_leaf;

  return y;
}


// ----------------------------------------------------
void ISOPRENE_CANOPY_FLUX( double *iso_pt)
{
  double iso_sun, iso_shade, iso_emission, zz;
  int j;

  iso_emission = 0.;

  for (j=1; j<=jtot; j++)
    {
      zz = (int)(delz * j);

      ISOPRENE_LEAF_FLUX(solar.quantum_sun[j], prof.sun_tleaf_filter[j], zz, &iso_sun);
      ISOPRENE_LEAF_FLUX(solar.quantum_shd[j], prof.shd_tleaf_filter[j], zz, &iso_shade);

      prof.iso_sun[j] = iso_sun;
      prof.iso_shd[j] = iso_shade;

      prof.sun_isopreneflux[j] = prof.dLAIdz[j] * (solar.prob_beam[j] * iso_sun);
      prof.shd_isopreneflux[j] = prof.dLAIdz[j] * (solar.prob_shd[j] * iso_shade);

      prof.isopreneflux[j]     = prof.sun_isopreneflux[j] + prof.shd_isopreneflux[j];

      iso_emission += prof.isopreneflux[j];
    }

  *iso_pt = iso_emission;

  return;
}


// ----------------------------------------------------
void ISOPRENE_LEAF_FLUX(double ppfd, double tlf, double zzz, double *isofluxpt)
{
  // leaf level algorithms of Guenther/Harley for oak leaves
  double t_lfk, light, c_t1, c_t2, isoref;

  t_lfk = tlf + TN0;

  // equation from original CANOAK
  light = 0.0027 * ppfd * 1.066 / sqrt(1. + 0.0027 * 0.0027 * ppfd * ppfd);
  c_t1 = exp(95000. * (t_lfk - 303.) / (rugc * t_lfk * 303.));
  c_t2 = 1. + exp(230000. * (303. - 314.) / (rugc * 303. * 314.));
  // nmol m-2 s-1
  isoref = 45.8 * (0.65 * zzz / ht + 0.35);

  /*
  // equation after adaptation from Jenifer Funk (Stanford)
  light = .0015 * ppfd * 1.1032 / sqrt(1. + .0015 * .0015 * ppfd * ppfd);
  c_t1 = exp(93960 * (t_lfk - 301.15) / (rugc * t_lfk * 301.15));
  c_t2 = 1 + exp(540600 * (t_lfk - 311.3) / (rugc * t_lfk * 301.15));
  // nmol m-2 s-1
  isoref = 32.1 * (.65 * zzz / ht + .35);
  */

  *isofluxpt = isoref * light * c_t1 / c_t2;

  return;
}


// ----------------------------------------------------
void LAI_TIME()
{
  // Evaluate how LAI and other canopy structural variables vary
  // with time

  long int I, J, II, JM1;

  double lai_z[sze];
  double TF,MU1,MU2,integr_beta;
  double dx,DX2,DX4,X,P_beta,Q_beta,F1,F2,F3;
  double beta_fnc[sze];
  double lai_freq_local[betasze];
  double ht_midpt_local[betasze];
  double cum_lai,sumlai,dff,XX;
  double cum_ht;
  double AA,DA,dff_Markov;
  double cos_AA,sin_AA,exp_diffuse;

  // seasonal update of model parameters

  if (time_var.days < time_var.leafout) // winter
    {
      time_var.lai = pai;
      // PAR wave band: after Norman (1979) and NASA report
      solar.par_reflect = par_reflect[1];
      solar.par_trans = par_trans[1];
      solar.par_soil_refl_dry = par_soil_refl_dry[1];
      solar.par_absorbed = (1. - solar.par_reflect - solar.par_trans);
      // NIR wave band: after Norman (1979) and NASA report
      solar.nir_reflect = nir_reflect[1];
      solar.nir_trans = nir_trans[1];
      solar.nir_soil_refl_dry = nir_soil_refl_dry[1];
      solar.nir_absorbed = (1. - solar.nir_reflect - solar.nir_trans);
    }

  if (time_var.days >= time_var.leafout && time_var.days < time_var.leaffull) // spring
    {
      time_var.lai = pai + (time_var.days - time_var.leafout) * (lai - pai) / (time_var.leaffull-time_var.leafout);
      // PAR wave band: after Norman (1979) and NASA report
      solar.par_reflect = par_reflect[2];
      solar.par_trans = par_trans[2];
      solar.par_soil_refl_dry = par_soil_refl_dry[2];
      solar.par_absorbed = (1. - solar.par_reflect - solar.par_trans);
      // NIR wave band: after Norman (1979) and NASA report
      solar.nir_reflect = nir_reflect[2];
      solar.nir_trans = nir_trans[2];
      solar.nir_soil_refl_dry = nir_soil_refl_dry[2];
      solar.nir_absorbed = (1. - solar.nir_reflect - solar.nir_trans);
    }

  if (time_var.days >= time_var.leaffull && time_var.days < time_var.leaffall) //summer
    {
      time_var.lai = lai;
      // PAR wave band: after Norman (1979) and NASA report
      solar.par_reflect = par_reflect[3];
      solar.par_trans = par_trans[3];
      solar.par_soil_refl_dry = par_soil_refl_dry[3];
      solar.par_absorbed = (1. - solar.par_reflect - solar.par_trans);
      // NIR wave band: after Norman (1979) and NASA report
      solar.nir_reflect = nir_reflect[3];
      solar.nir_trans = nir_trans[3];
      solar.nir_soil_refl_dry = nir_soil_refl_dry[3];
      solar.nir_absorbed = (1. - solar.nir_reflect - solar.nir_trans);
    }

  if (time_var.days >= time_var.leaffall && time_var.days < time_var.leaffallcomplete) //fall
    {
      time_var.lai = lai - (time_var.days - time_var.leaffall) * (lai-pai) / (time_var.leaffallcomplete-time_var.leaffall); //changed
      // PAR wave band: after Norman (1979) and NASA report
      solar.par_reflect = par_reflect[4];
      solar.par_trans = par_trans[4];
      solar.par_soil_refl_dry = par_soil_refl_dry[4];
      solar.par_absorbed = (1. - solar.par_reflect - solar.par_trans);
      // NIR wave band: after Norman (1979) and NASA report
      solar.nir_reflect = nir_reflect[4];
      solar.nir_trans = nir_trans[4];
      solar.nir_soil_refl_dry = nir_soil_refl_dry[4];
      solar.nir_absorbed = (1. - solar.nir_reflect - solar.nir_trans);
    }

  if (time_var.days >= time_var.leaffallcomplete) //winter
    {
      time_var.lai = pai;
      // PAR wave band: after Norman (1979) and NASA report
      solar.par_reflect = par_reflect[1];
      solar.par_trans = par_trans[1];
      solar.par_soil_refl_dry = par_soil_refl_dry[1];
      solar.par_absorbed = (1. - solar.par_reflect - solar.par_trans);
      // NIR wave band: after Norman (1979) and NASA report
      solar.nir_reflect = nir_reflect[1];
      solar.nir_trans = nir_trans[1];
      solar.nir_soil_refl_dry = nir_soil_refl_dry[1];
      solar.nir_absorbed = (1. - solar.nir_reflect - solar.nir_trans);
    }

  // for Nate McDowell's juniper site read LAI instead of diffuse PAR
  if (extra_nate == 1)
    time_var.lai = __max(input.lai_up+input.lai_down, pai);
  // height of mid point of layer scaled to height of forest

  for (I=1; I<=betasze-1; I++)
    ht_midpt_local[I] = ht_midpt[I];

  if (time_var.days >= time_var.leafout && time_var.days <= time_var.leaffallcomplete) // use LAI distribution
    {
      // lai of the layers at the midpoint of height
      for (I=1; I<=betasze-1; I++)
	lai_freq_local[I] = lai_freq[I]*lai;
    }
  else // use PAI distribution
    {
      for (I=1; I<=betasze-1; I++)
	lai_freq_local[I] = pai_freq[I]*pai;
    }

  /*
    Beta distribution
    f(x) = x^(p-1) (1-x)^(q-1) / B(v,w)
    B(v,w) = int from 0 to 1 x^(p-1) (1-x)^(q-1) dx
    p = mean{[mean(1-mean)/var]-1}
    q =(1-mean){[mean(1-mean)/var]-1}
  */
  TF = 0.;
  MU1 = 0.;
  MU2 = 0.;
  integr_beta = 0.;
  /*
    ENTER THE HEIGHT AT THE MIDPOINT OF THE LAYER
  */
  for (I=1; I<=betasze-1; I++)
    {
      // Normalize height
      ht_midpt_local[I] /= ht;
      // TOTAL F IN EACH LAYER, TF SHOULD SUM TO EQUAL lai
      TF += lai_freq_local[I];
      // weighted mean lai
      MU1 += (ht_midpt_local[I] * lai_freq_local[I]);
      // weighted variance
      MU2 += (ht_midpt_local[I] * ht_midpt_local[I] * lai_freq_local[I]);
    } // next I

  // normalize mu by lai
  MU1 /= TF;
  MU2 /= TF;

  // compute Beta parameters
  P_beta = MU1 * (MU1 - MU2) / (MU2 - MU1 * MU1);
  Q_beta = (1. - MU1) * (MU1 - MU2) / (MU2 - MU1 * MU1);
  P_beta -= 1.;
  Q_beta -= 1.;

  /*
    ' integrate Beta function, with Simpson's Approx.
    '
    ' The boundary conditions are level 1 is height of ground
    ' and level jtot+1 is height of canopy. Layer 1 is between
    ' height levels 1 and 2. Layer jtot is between levels
    ' jtot and jtot+1
    '
    ' Thickness of layer
  */

  dx = 1. / jtot;

  DX2 = dx / 2.;
  DX4 = dx / 4.;
  X = DX4;

  F2 = pow(X, P_beta) * pow((1.-X), Q_beta);
  X += DX4;
  F3 = pow(X, P_beta) *pow((1.-X), Q_beta);

  // start integration at lowest boundary


  beta_fnc[1] = DX4 * (4. * F2 + F3) / 3.;
  integr_beta += beta_fnc[1];

  JM1=jtot-1;

  for (I = 2; I<=JM1; I++)
    {
      F1 = F3;
      X += DX2;
      F2 = pow(X, P_beta) * pow((1. - X), Q_beta);
      X += DX2;
      F3 = pow(X, P_beta) * pow((1. - X), Q_beta);
      beta_fnc[I] = DX2 * (F1 + 4. * F2 + F3) / 3.;
      integr_beta += beta_fnc[I];
    }

  F1 = F3;
  X += DX4;
  F2 = pow(X, P_beta) * pow((1. - X),Q_beta);

  // compute integrand at highest boundary

  beta_fnc[jtot] = DX4 * (F1 + 4. * F2) / 3.;
  integr_beta += beta_fnc[jtot];
  /*
    ' lai_z IS THE LEAF AREA AS A FUNCTION OF Z
    '
    ' beta_fnc is the pdf for the interval dx
  */

  lai_z[1] = beta_fnc[1] * time_var.lai / integr_beta;

  for (I=2; I<=JM1; I++)
    lai_z[I] = beta_fnc[I] * time_var.lai / integr_beta;

  lai_z[jtot] = beta_fnc[jtot] * time_var.lai / integr_beta;

  cum_ht = 0;
  cum_lai = 0;


  for (I=1; I<=jtot; I++)
    {
      /*
        ' re-index layers of lai_z.
        ' layer 1 is between ground and 1st level
        ' layer jtot is between level jtot and top of canopy (jtot+1)
      */
      cum_ht += delz;
      cum_lai += lai_z[I];
      // use prof.dLAIdz for radiative transfer model
      prof.dLAIdz[I] = lai_z[I];
    } // next I

  // for Nate McDowell's juniper site, hardcode lai distribution for 40 layers
  if (extra_nate == 1)
    {
      // equall distribution in lowest 50cm = 8 layers -> nup=9
      for (I=1; I<=nup-1; I++)
	prof.dLAIdz[I] = __max(input.lai_down/((double) (nup-1)), pai/jtot);
      // pear shaped profile of upper canopy
      // start tapering from ca. 2m = layer 31
      for (I=nup; I<=jtot-nup-1; I++)
	prof.dLAIdz[I]  = 1.0;
      for (I=jtot-nup; I<=jtot; I++)
	//MC20190615 prof.dLAIdz[I] = 1.0 - (1./(((float) nup)+1.))*(nup-(jtot-I));
	prof.dLAIdz[I] = 1.0 - (1./(((double) nup)+1.))*(nup-(jtot-I));
      // normalise pear shape
      sumlai = 0.;
      for (J=nup; J<=jtot; J++)
	sumlai += prof.dLAIdz[J];
      for (J=nup; J<=jtot; J++)
	prof.dLAIdz[J] = __max(prof.dLAIdz[J]/sumlai * input.lai_up, pai/jtot);
      // sync time_var.lai with dLAIdz
      time_var.lai = 0.;
      for (J=1; J<=jtot; J++)
	time_var.lai += prof.dLAIdz[J];
    }
  else
    {
      // normalise total LAI including layers with PAI
      sumlai = 0;
      for (J=1; J<=jtot; J++)
	sumlai += prof.dLAIdz[J];
      for (J=1; J<=jtot; J++)
	prof.dLAIdz[J] = prof.dLAIdz[J]*time_var.lai/sumlai;
    }

  // PAI
  for (I=1; I<=jtot; I++)
    prof.dPAIdz[I] = prof.dLAIdz[I]/time_var.lai*pai;

  G_FUNC_DIFFUSE(); // Direction cosine for the normal between the mean

  for (J=1; J<=jtot; J++)
    {
      /*
        '
        ' compute the probability of diffuse radiation penetration through the
        ' hemisphere. This computation is not affected by penubra
        ' since we are dealing only with diffuse radiation from a sky
        ' sector.
        '
        ' The probability of beam penetration is computed with a
        ' Markov distribution.
      */
      dff = prof.dLAIdz[J];//+ prof.dPAIdz[J];
      XX = 0;
      AA = .087;
      DA = .1745;
      // The leaf clumping coefficient. From Chason et al. 1990 and studies on WBW
      dff_Markov = dff*markov;
      for (II=1; II<=9;II++)
        {
          cos_AA = cos(AA);
          sin_AA = sin(AA);

          // probability of photon transfer through a sky section
          // for clumped foliage and Markov model
          // for spherical distribution
          // exp_diffuse = exp(-DFF * prof.Gfunc_sky(J, II) / cos_AA)
          exp_diffuse = exp(-dff_Markov * prof.Gfunc_sky[J][II] / cos_AA);

          XX += (cos_AA * sin_AA * exp_diffuse);
          AA += DA;
        } // next II

      /*
        'Itegrated probability of diffuse sky radiation penetration
        'for each layer
      */
      solar.exxpdir[J] = 2. * XX * DA;
      if (solar.exxpdir[J] > 1.)
	solar.exxpdir[J] = 1.;
    } // next J

  return;
}


// ----------------------------------------------------
void NIR()
{

  /*
    ------------------------------------------------------------
    SUBROUTINE NIR

    This subroutine computes the flux density of direct and diffuse
    radiation in the near infrared waveband. The Markov model is used
    to compute the probability of beam penetration through clumped foliage.

    The algorithms of Norman (1979) are used.

    Norman, J.M. 1979. Modeling the complete crop canopy.
    Modification of the Aerial Environment of Crops.
    B. Barfield and J. Gerber, Eds. American Society of Agricultural Engineers, 249-280.
    -----------------------------------------------------------------
  */

  long int J, JM1, JJ, JJP1, ITER,IREP;

  double fraction_beam;
  double SUP[sze], SDN[sze], transmission_layer[sze], reflectance_layer[sze], beam[sze];
  double TBEAM[sze];
  double ADUM[sze];
  double exp_direct,sumlai, dff;
  double TLAY2,nir_normal, NSUNEN;
  double llai,NIRTT,DOWN, UP;

  solar.nir_total = solar.nir_beam + solar.nir_diffuse;

  fraction_beam = solar.nir_beam / (solar.nir_beam + solar.nir_diffuse);
  beam[jtot1] = fraction_beam;
  TBEAM[jtot1] = fraction_beam;

  solar.nir_dn[jtot1] = 1. - fraction_beam;

  if (solar.nir_total <= 1. || solar.sine_beta <= 0.05)
    goto NIRNIGHT;

  SDN[1] = 0;

  /*

  Compute probability of penetration for direct and
  diffuse radiation for each layer in the canopy

  Level 1 is the soil surface and level jtot1 is the
  top of the canopy. layer 1 is the layer above
  the soil and layer jtot is the top layer.
  */

  sumlai = 0.;

  for (J=2; J<=jtot1;J++)
    {
      JJ = jtot1 - J + 1;

      /*
        diffuse NIR reflected by each layer

        compute the probability of diffuse radiation penetration through the
        hemisphere. this computation is not affected by penubra
        since we are dealing only with diffuse radiation from a sky
        sector.

        The probability of beam penetration is computed with a
        negative binomial distriubution for LAI.

        Radiation attenuation is a function of leaf and woody
        biomass
      */

      sumlai += prof.dLAIdz[JJ]; // + prof.dPAIdz[JJ]

      /*
        'Itegrated probability of diffuse sky radiation penetration
        'for each layer
        '
        EXPDIF[JJ] is computed in PAR and can be used in NIR and IRFLUX
      */

      reflectance_layer[JJ] = (1. - solar.exxpdir[JJ]) * solar.nir_reflect;


      // DIFFUSE RADIATION TRANSMITTED THROUGH LAYER

      transmission_layer[JJ] = (1. - solar.exxpdir[JJ]) * solar.nir_trans + solar.exxpdir[JJ];
    } // next J


  // COMPUTE THE PROBABILITY OF beam PENETRATION

  for (J=2; J<=jtot1; J++)
    {
      JJ = jtot1 - J + 1;
      JJP1 = JJ + 1;

      // Probability of beam penetration.


      dff = prof.dLAIdz[JJ]; /* '+ prof.dPAIdz[JJ] */

      exp_direct = exp(-dff * markov*prof.Gfunc_solar[JJ] / solar.sine_beta);

      // PEN1 = exp(-llai * prof.Gfunc_solar[JJ] / solar.sine_beta)
      // exp_direct = exp(-DFF * prof.Gfunc_solar[JJ] / solar.sine_beta)

      // Beam transmission

      beam[JJ] = beam[JJP1] * exp_direct;

      TBEAM[JJ] = beam[JJ];


      SUP[JJP1] = (TBEAM[JJP1] - TBEAM[JJ]) * solar.nir_reflect;

      SDN[JJ] = (TBEAM[JJP1] - TBEAM[JJ]) * solar.nir_trans;
    } // next J

  /*
    initiate scattering using the technique of NORMAN (1979).
    scattering is computed using an iterative technique
  */

  SUP[1] = TBEAM[1] * solar.nir_soil_refl;
  ADUM[1] = solar.nir_soil_refl;

  for (J = 2; J<=jtot1; J++)
    {
      JM1 = J - 1;
      TLAY2 = transmission_layer[JM1] * transmission_layer[JM1];
      ADUM[J] = ADUM[JM1] * TLAY2 / (1. - ADUM[JM1] * reflectance_layer[JM1]) + reflectance_layer[JM1];
    } // NEXT J

  for (J = 1; J<=jtot; J++)
    {
      JJ = jtot - J + 1;
      JJP1 = JJ + 1;

      solar.nir_dn[JJ] = solar.nir_dn[JJP1] * transmission_layer[JJ] / (1. - ADUM[JJP1] * reflectance_layer[JJ]) + SDN[JJ];
      solar.nir_up[JJP1] = ADUM[JJP1] * solar.nir_dn[JJP1] + SUP[JJP1];
    }

  // lower boundary: upward radiation from soil


  solar.nir_up[1] = solar.nir_soil_refl * solar.nir_dn[1] + SUP[1];

  /*
    Iterative calculation of upward diffuse and downward beam +
    diffuse NIR to compute scattering
  */

  ITER = 0;
  IREP = 1;

  ITER += 1;

  while (IREP == 1)
    {

      IREP=0;
      for (J = 2; J<=jtot1;J++)
        {
          JJ = jtot1 - J + 1;
          JJP1 = JJ + 1;
          DOWN = transmission_layer[JJ] * solar.nir_dn[JJP1] + solar.nir_up[JJ] * reflectance_layer[JJ] + SDN[JJ];

          if ((fabs(DOWN - solar.nir_dn[JJ])) > .01) IREP = 1;

          solar.nir_dn[JJ] = DOWN;
        }

      // upward radiation at soil is reflected beam and downward diffuse

      solar.nir_up[1] = (solar.nir_dn[1] + TBEAM[1]) * solar.nir_soil_refl;

      for (JJ = 2; JJ <=jtot1;JJ++)
        {
          JM1 = JJ - 1;

          UP = reflectance_layer[JM1] * solar.nir_dn[JJ] + solar.nir_up[JM1] * transmission_layer[JM1] + SUP[JJ];

          if ((fabs(UP - solar.nir_up[JJ])) > .01) IREP = 1;

          solar.nir_up[JJ] = UP;
        }

    }


  // Compute NIR flux densities

  llai = sumlai;

  for (J = 1; J<=jtot1;J++)
    {
      llai -= prof.dLAIdz[J]; // decrement LAI

      // upward diffuse NIR flux density, on the horizontal
      solar.nir_up[J] *= solar.nir_total;

      if (solar.nir_up[J] < 0.1)
	solar.nir_up[J] = 0.1;

      // downward beam NIR flux density, incident on the horizontal
      solar.beam_flux_nir[J] = beam[J] * solar.nir_total;

      if (solar.beam_flux_nir[J] < 0.1)
	solar.beam_flux_nir[J] = 0.1;

      // downward diffuse radiaiton flux density on the horizontal
      solar.nir_dn[J] *= solar.nir_total;

      if (solar.nir_dn[J] < 0.1)
	solar.nir_dn[J] = 0.1;

      // total downward NIR, incident on the horizontal
      NIRTT = solar.beam_flux_nir[J] + solar.nir_dn[J];
    } // next J

  for (J = 1; J<=jtot; J++)
    {
      // normal radiation on sunlit leaves
      if (solar.sine_beta > 0.05)
        nir_normal = solar.nir_beam * prof.Gfunc_solar[J] / solar.sine_beta;
      else
        nir_normal = 0.;

      NSUNEN = nir_normal * solar.nir_absorbed;

      /*
        Diffuse radiation received on top and bottom of leaves
        drive photosynthesis and energy exchanges
      */
      solar.nir_shd[J] = (solar.nir_dn[J] + solar.nir_up[J]);

      // absorbed radiation, shaded
      solar.nir_shd[J] *= solar.nir_absorbed;

      // plus diffuse component
      solar.nir_sun[J] = NSUNEN + solar.nir_shd[J];
    } // next J

 NIRNIGHT: // jump to here at night since fluxes are zero

  if (solar.nir_total <= 1. || solar.sine_beta <= 0.05)
    {
      for (J = 1; J<=jtot; J++)
        {
          solar.nir_up[J] = 0.;
          solar.nir_dn[J] = 0.;
          solar.nir_shd[J] = 0.;
          solar.nir_sun[J] = 0.;
          solar.beam_flux_nir[J] = 0.;
        }
      solar.nir_total = 0.;
      solar.nir_up[jtot1] = 0.;
      solar.nir_dn[jtot1] = 0.;
      solar.beam_flux_nir[jtot1] = 0.;
    }

  return;
}


// ----------------------------------------------------
void PAR()
{
  /*
    -------------------------------------------------------------
    SUBROUTINE PAR

    This subroutine computes the flux densities of direct and
    diffuse radiation using the measured leaf distrib.
    We apply the Markov model to account for clumping of leaves.

    The model is based on the scheme of Norman.

    Norman, J.M. 1979. Modeling the complete crop canopy.
    Modification of the Aerial Environment of Crops.
    B. Barfield and J. Gerber, Eds. American Society of Agricultural Engineers, 249-280.

    -----------------------------------------------------------

  */

  double SUP[sze], SDN[sze], transmission_layer[sze], reflectance_layer[sze], beam[sze];
  double ADUM[sze], TBEAM[sze];
  double fraction_beam,sumlai;
  double PEN1, PEN2, dff, par_normal_quanta;
  double exp_direct,QU,TLAY2;
  double par_normal_abs_energy, par_normal_abs_quanta;
  double DOWN, UP;

  long int J, JJP1;
  long int I, JM1, JJ, par_total;
  long int IREP, ITER;


  if (solar.sine_beta <= 0.05) goto
    par_night;

  fraction_beam = solar.par_beam / input.parin;

  beam[jtot1] = fraction_beam;
  TBEAM[jtot1] = fraction_beam;

  SDN[1] = 0;

  /*
    Compute probability of penetration for direct and
    diffuse radiation for each layer in the canopy

    Level 1 is the soil surface and level jtot1 is the
    top of the canopy. layer 1 is the layer above
    the soil and layer jtot is the top layer.
  */
  PEN1 = 1.;

  for (J=1; J<=jtot;J++)
    {
      JJ = jtot1 - J;

      /*
        Diffuse PAR reflected by each layer

        compute the probability of diffuse radiation penetration through the
        hemisphere. this computation is not affected by penubra
        since we are dealing only with diffuse radiation from a sky
        sector.

        The probability of beam penetration is computed with a
        negative binomial distriubution.
      */

      dff = prof.dLAIdz[JJ]; /* + prof.dPAIdz[JJ] */

      reflectance_layer[JJ] = (1. - solar.exxpdir[JJ]) * solar.par_reflect;

      // DIFFUSE RADIATION TRANSMITTED THROUGH LAYER

      transmission_layer[JJ] = (1. - solar.exxpdir[JJ]) * solar.par_trans + solar.exxpdir[JJ];

    } // next J

  /*
    'COMPUTE THE PROBABILITY OF beam PENETRATION
  */
  sumlai = 0.;
  for (J=1; J<=jtot;J++)
    {
      JJ = jtot1 - J;
      JJP1 = JJ + 1;

      /*
        Probability of beam penetration. This is that which
        is not umbral. Hence, this radiation
        is attenuated by the augmented leaf area: DF+PA.
      */
      dff = prof.dLAIdz[JJ]; /* + prof.dPAIdz[JJ] */

      sumlai += dff;

      exp_direct = exp(-dff*markov*prof.Gfunc_solar[JJ]/ solar.sine_beta);

      PEN2 = exp(-sumlai*markov*prof.Gfunc_solar[JJ]/ solar.sine_beta);


      /* lai Sunlit and shaded */

      prof.sun_lai[JJ] = solar.sine_beta * (1 - PEN2)/ (markov*prof.Gfunc_solar[JJ]);

      prof.shd_lai[JJ] = sumlai - prof.sun_lai[JJ];

      /* note that the integration of the source term time solar.prob_beam with respect to
         leaf area will yield the sunlit leaf area, and with respect to solar.prob_shd the
         shaded leaf area.


         In terms of evaluating fluxes for each layer

         Fcanopy = sum {fsun psun + fshade pshade} (see Leuning et al. Spitters et al.)

         psun is equal to exp(-lai G markov/sinbet)

         pshade = 1 - psun


      */

      solar.prob_beam[JJ] = markov*PEN2;

      if (solar.prob_beam[JJ] <= 0.) PEN1 = 0;


      // probability of beam

      beam[JJ] = beam[JJP1] * exp_direct;

      QU = 1.0 - solar.prob_beam[JJ];

      if (QU > 1.) QU = 1.;

      if (QU < 0.) QU = 0.;


      // probability of umbra

      solar.prob_shd[JJ] = QU;

      TBEAM[JJ] = beam[JJ];


      // beam PAR that is reflected upward by a layer

      SUP[JJP1] = (TBEAM[JJP1] - TBEAM[JJ]) * solar.par_reflect;


      // beam PAR that is transmitted downward


      SDN[JJ] = (TBEAM[JJP1] - TBEAM[JJ]) * solar.par_trans;

    } // next J
  /*
    initiate scattering using the technique of NORMAN (1979).
    scattering is computed using an iterative technique.

    Here Adum is the ratio up/down diffuse radiation.

  */
  SUP[1] = TBEAM[1] * solar.par_soil_refl;

  solar.par_down[jtot1] = 1.0 - fraction_beam;
  ADUM[1] = solar.par_soil_refl;
  for (J = 2,JM1=1;J<=jtot1;J++,JM1++)
    {
      TLAY2 = transmission_layer[JM1] * transmission_layer[JM1];
      ADUM[J] = ADUM[JM1] * TLAY2 / (1. - ADUM[JM1] * reflectance_layer[JM1]) + reflectance_layer[JM1];
    } /* NEXT J */

  for (J=1; J<=jtot; J++)
    {
      JJ = jtot - J + 1;
      JJP1 = JJ + 1;
      solar.par_down[JJ] = solar.par_down[JJP1] * transmission_layer[JJ] / (1. - ADUM[JJP1] * reflectance_layer[JJ]) + SDN[JJ];
      solar.par_up[JJP1] = ADUM[JJP1] * solar.par_down[JJP1] + SUP[JJP1];
    } // next J

  // lower boundary: upward radiation from soil

  solar.par_up[1] = solar.par_soil_refl * solar.par_down[1] + SUP[1];

  /*
    Iterative calculation of upward diffuse and downward beam +
    diffuse PAR.

    This section has been commented out for the negative binomial
    model. It seems not to apply and incorrectly calculates
    scattering. When I ignore this section, I get perfect
    agreement between measured and calculated Rn.

  */

  // Scattering

  ITER = 0;
  IREP=1;

  while (IREP == 1)
    {
      IREP = 0;

      ITER += 1;

      for (J=2; J<=jtot1; J++)
        {
          JJ = jtot1 - J + 1;
          JJP1 = JJ + 1;
          DOWN = transmission_layer[JJ] * solar.par_down[JJP1] + solar.par_up[JJ] * reflectance_layer[JJ] + SDN[JJ];

          if ((fabs(DOWN - solar.par_down[JJ])) > .01) IREP = 1;

          solar.par_down[JJ] = DOWN;
        } // next J

      // upward radiation at soil is reflected beam and downward diffuse */

      solar.par_up[1] = (solar.par_down[1] + TBEAM[1]) * solar.par_soil_refl;

      for (JJ = 2; JJ <=jtot1; JJ++)
        {
          JM1 = JJ - 1;
          UP = reflectance_layer[JM1] * solar.par_down[JJ] + solar.par_up[JM1] * transmission_layer[JM1] + SUP[JJ];

          if ((fabs(UP - solar.par_up[JJ])) > .01) IREP = 1;

          solar.par_up[JJ] = UP;
        } // next JJ

    }

  // Compute flux density of PAR

  for (J = 1;J<=jtot1;J++)
    {
      // upward diffuse PAR flux density, on the horizontal
      solar.par_up[J] *= input.parin;

      if (solar.par_up[J] < 0.001)
	solar.par_up[J] = 0.001;

      // downward beam PAR flux density, incident on the horizontal
      solar.beam_flux_par[J] = beam[J] * input.parin;

      if (solar.beam_flux_par[J] < 0.001) 
	solar.beam_flux_par[J] = 0.001;

      // Downward diffuse radiatIon flux density on the horizontal
      solar.par_down[J] *= input.parin;

      if (solar.par_down[J] < 0.001)
	solar.par_down[J] = 0.001;

      // Total downward PAR, incident on the horizontal
      par_total = (int) solar.beam_flux_par[J] + (int) solar.par_down[J];
    } // next J

  if (solar.par_beam < 0.001)
    solar.par_beam = 0.001;

  // PSUN is the radiation incident on the mean leaf normal

  for (JJ = 1;JJ<=jtot;JJ++)
    {
      if (solar.sine_beta > 0.05)
        par_normal_quanta = solar.par_beam * prof.Gfunc_solar[JJ] / (solar.sine_beta); // PAR received normal to a leaf on a sunlit spot
      else
        par_normal_quanta = 0.;

      // amount of energy absorbed by sunlit leaf */

      par_normal_abs_energy = par_normal_quanta * solar.par_absorbed / 4.6; // W m-2
      par_normal_abs_quanta = par_normal_quanta * solar.par_absorbed; // umol m-2 s-1

      /*
        Convert PAR to W m-2 for energy balance computations
        but remember umol m-2 s-1 is needed for photosynthesis.
        Average fluxes for the values on the top and bottom of
        each layer that drives energy fluxes.

        Energy balance computations are on the basis of
        absorbed energy. Harley's photosynthesis
        parameterizations are on the basis of incident
        PAR
      */
      solar.quantum_shd[JJ] = (solar.par_down[JJ] + solar.par_up[JJ])*solar.par_absorbed; /* umol m-2 s-1 */
      solar.quantum_sun[JJ] = solar.quantum_shd[JJ] + par_normal_abs_quanta;


      // calculate absorbed par
      solar.par_shd[JJ] = solar.quantum_shd[JJ] / 4.6; /* W m-2 */

      /*
        solar.par_sun is the total absorbed radiation on a sunlit leaf,
        which consists of direct and diffuse radiation
      */
      solar.par_sun[JJ] = par_normal_abs_energy + solar.par_shd[JJ];
    } // next JJ

 par_night:

  if (solar.sine_beta <= 0.05)
    {
      for (I = 1; I<=jtot; I++)
        {
          solar.prob_shd[I] = 1.;
          solar.prob_beam[I] = 0.;
          solar.par_up[I] = 0.;
          solar.par_down[I] = 0.;
          solar.par_sun[I] = 0.;
          solar.par_shd[I] = 0.;
          solar.beam_flux_par[I] = 0.;
          solar.quantum_shd[I] = 0.;
          solar.quantum_sun[I] = 0.;
          prof.sun_lai[I] = 0.;
          prof.shd_lai[I] = prof.dLAIdz[I];
        }
      solar.par_beam = 0.;
      solar.par_up[jtot1] = 0.;
      solar.par_down[jtot1] = 0.;
      solar.beam_flux_par[jtot1] = 0.;
    }
  return;
}


// ----------------------------------------------------
void FRICTION_VELOCITY()
{
  // this subroutine updates ustar and stability corrections
  // based on the most recent H and z/L values

  double xzl, logprod, phim;

  met.zl = -(vonKarman*Gravity*met.H_filter*(zm-zd))
    /(met.air_density*1005.*met.T_Kelvin
      *met.ustar_filter*met.ustar_filter*met.ustar_filter); // z/L

  // restrict z/L within reasonable bounds to minimize numerical errors
  if (met.zl > 0.25)
    met.zl = 0.25;
  if (met.zl < -3.)
    met.zl = -3.;

  // calculation based on Businger et a. 1977 and Hoegstroem 1988 (see Foken 2003)
  if (met.zl < 0.)
    {
      xzl = pow(1.-met.zl*16., 0.25);
      logprod = (1.+xzl*xzl)/2. * (1.+xzl)/2.*(1.+xzl)/2.;
      if (logprod > 0. || logprod < 100000.)
        logprod = log(logprod);
      else
        logprod = -1.;
      phim = logprod - 2.*atan(xzl) + PI/2.;
    }
  else
    phim = -5.*met.zl;

  if ((log((zm-zd)/z0)-phim) > 0.)
    met.ustar = vonKarman*input.wnd/(log((zm-zd)/z0)-phim);
  else
    met.ustar = met.ustar;

  if (met.ustar > 2.)
    met.ustar = vonKarman*input.wnd/log((zm-zd)/z0);

  if (met.ustar < 0.02)
    met.ustar = 0.02;

  return;
}


// ----------------------------------------------------
void BOUNDARY_RESISTANCE(double zzz, double TLF, double cws, int JLAY)
{
  /*
    BOUNDARY_RESISTANCE

    This subroutine computes the leaf boundary layer
    resistances for heat, vapor and CO2 (s/m).

    Flat plate theory is used, as discussed in Schuepp (1993) and
    Grace and Wilson (1981).

    We consider the effects of turbulent boundary layers and sheltering.
    Schuepp's review shows a beta factor multiplier is necessary for SH in
    flows with high turbulence. The concepts and theories used have been
    validated with our work on HNO3 transfer to the forest.

    Schuepp. 1993 New Phytologist 125: 477-507

    Diffusivities have been corrected using the temperature/Pressure algorithm in Massman (1998)
  */

  double Re, Re5, Re8; // Reynolds numbers
  double Sh_heat, Sh_vapor, Sh_CO2; // Sherwood numbers
  double graf, GR25; // Grasshof numbers
  double deltlf;
  double Res_factor;//, ddia, dRe;
  double nnu_T_P, ddh_T_P, ddv_T_P, ddc_T_P, T_kelvin;

  /* 'Difference between leaf and air temperature */

  deltlf = (TLF - prof.tair_filter[JLAY]);

  T_kelvin = prof.tair_filter[JLAY] + TN0;

  if (deltlf > 0.)
    graf = non_dim.grasshof * deltlf / T_kelvin;
  else
    graf = 0.;

  nnu_T_P = nnu * (1013./input.press_mb) * pow((T_kelvin/TN0), 1.81);

  prof.u[JLAY] = UZ(zzz);

  Re = lleaf * prof.u[JLAY] / nnu_T_P;

  if (Re > 0.)
    Re5 = sqrt(Re);
  else
    {
      puts("bad RE in RESHEAT\n");
      Re5=100.;
    }

  Re8 = pow(Re, 0.8);

  if (Re > 14000.)
    {
      Res_factor = 0.036*Re8*betfact;
      /*
        turbulent boundary layer
        SH = .036 * Re8 * pr33*betfact;
        SHV = .036 * Re8 * sc33*betfact;
        SHCO2 = .036 * Re8 * scc33*betfact;
      */
      Sh_heat = Res_factor * non_dim.pr33;
      Sh_vapor = Res_factor * non_dim.sc33;
      Sh_CO2 = Res_factor * non_dim.scc33;
    }
  else
    {
      Res_factor = 0.66*Re5*betfact;
      /*
        laminar sublayer
        SH = .66 * Re5 * pr33*betfact;
        SHV = .66 * Re5 * sc33*betfact;
        SHCO2 = .66 * Re5 * scc33*betfact;
      */
      Sh_heat = Res_factor * non_dim.pr33;
      Sh_vapor = Res_factor * non_dim.sc33;
      Sh_CO2 = Res_factor * non_dim.scc33;
      if (cws > 0.)
	Sh_vapor = 0.66 * non_dim.sc33 * pow(Re, 0.4);
      //   Sh_vapor = 0.66 * non_dim.sc33 * pow(Re, 0.4) * betfact;
    }


  // If there is free convection

  if (graf/(Re*Re) > 1.)
    {
      // Compute Grashof number for free convection
      if (graf < 100000.)
        GR25 = 0.5 * pow(graf, 0.25);
      else
        GR25 = 0.13 * pow(graf, 0.33);
      Sh_heat = non_dim.pr33 * GR25;
      Sh_vapor = non_dim.sc33 * GR25;
      Sh_CO2 = non_dim.scc33 * GR25;
    }

  // Correct diffusivities for temperature and pressure

  ddh_T_P = ddh * (1013./input.press_mb) * pow((T_kelvin/TN0), 1.81);
  ddv_T_P = ddv * (1013./input.press_mb) * pow((T_kelvin/TN0), 1.81);
  ddc_T_P = ddc * (1013./input.press_mb) * pow((T_kelvin/TN0), 1.81);

  bound_layer_res.heat = lleaf/(ddh_T_P * Sh_heat);
  bound_layer_res.vapor = lleaf/(ddv_T_P * Sh_vapor);
  bound_layer_res.co2 = lleaf / (ddc_T_P * Sh_CO2);

  return;
}


// ----------------------------------------------------
void RNET()
{
  double ir_shade, q_sun, q_shd;
  int j;
  /*
    Radiation layers go from jtot+1 to 1. jtot+1 is the top of the
    canopy and level 1 is at the soil surface.

    Energy balance and photosynthesis are performed for vegetation
    between levels and based on the energy incident to that level
  */
  for (j=1; j<=jtot; j++)
    {
      // Infrared radiation on leaves
      ir_shade = solar.ir_dn[j] + solar.ir_up[j];
      ir_shade = ir_shade * ep;
      /*
        Available energy on leaves for evaporation.
        Values are average of top and bottom levels of a layer.
        The index refers to the layer. So layer 3 is the average
        of fluxes at level 3 and 4. Level 1 is soil and level
        j+1 is the top of the canopy. layer jtot is the top layer

        Sunlit, shaded values
      */
      q_sun = solar.par_sun[j] + solar.nir_sun[j] + ir_shade;
      q_shd = solar.par_shd[j] + solar.nir_shd[j] + ir_shade;
      solar.rnet_sun[j] = q_sun;
      solar.rnet_shd[j] = q_shd;
    } // next j

  return;
}


// ----------------------------------------------------
void SET_SOIL_LITTER_CAPACITY()
{
  // Compute layer heat capacities and conductivities
  // as a function of texture, bulk density and soil moisture
  // Reference: Campbell, Soil physics with Basic (1985, Chapter 4)
  int i;
  double rho_0 = 1300.; // wet density of std organic matter [kg m-3]
  double c_0 = 2511600.; // heat capacity of std organic matter [J m-3 K-1]
  double c_w = 4185000.; // heat capacity of water [J m-3 K-1]
  double c_s = 2000000.; // heat capacity of soil matter [J m-3 K-1]
  // soil
  double C1, C2, C3, C4, C5;
  double rho_local;

  // Soil
  for (i=1; i<=soil.n_soil-1; i++) // as in Campbell (1985)
    {
      rho_local = (1.-soil.gravel[i]/100.)*soil.rho[i] + soil.gravel[i]/100.*2.65;

      soil.cp_soil[i] =
        ((c_s*rho_local/2.65) + c_w*soil.theta[i][1]*(1.-soil.gravel[i]/100.))
        * (soil.z_soil[i+1] - soil.z_soil[i-1]) / (2.*soil.temperature_dt);

      C1 = 0.65 - 0.78 * rho_local + 0.6 * rho_local*rho_local; //old=0.6 instead of 0.9
      C2 = 1.06 * rho_local; //old =1.06
      C3 = 1.00 + 2.6 / sqrt(soil.clay[i]);
      C4 = 0.03 + 0.1 * rho_local*rho_local;
      C5 = 4.;
      soil.k_conductivity_soil[i] = (C1 + 
					    C2 * soil.theta[i][1]*(1.-soil.gravel[i]/100.) - 
					    (C1 - C4) * exp(-pow(C3*soil.theta[i][1]*(1.-soil.gravel[i]/100.),C5)))
        / (soil.z_soil[i+1] - soil.z_soil[i]);
      /*// ISOLSM formulation
	C1 = (8.80*soil.sand[i]+2.92*soil.clay[i])/(soil.sand[i]+soil.clay[i]);
	C2 = 0.15;
	C3 = (pow(C1,1.-soil.theta_s[i]) * pow(0.6,soil.theta[i][1]) - C2 )
	* soil.theta[i][1]/soil.theta_s[i] + C2;
	soil.k_conductivity_soil[i] = C3 / (soil.z_soil[i+1] - soil.z_soil[i]);*/
    }
    /* printf("LC01 %20.14f %20.14f\n", rho_local[1], rho_local[soil.n_soil-1]); */

  // Last layer is different (not in Campbell (1985) who has an extra layer underneath the last soil layer)
  // Take thickness of last soil layer instead of mean thickness of two layers
  i = soil.n_soil;
  rho_local = (1.-soil.gravel[i]/100.)*soil.rho[i] + soil.gravel[i]/100.*2.65;

  soil.cp_soil[i] =
    ((c_s*rho_local/2.65) + c_w*soil.theta[i][1]*(1.-soil.gravel[i]/100.))
    * (soil.z_soil[i] - soil.z_soil[i-1]) / soil.temperature_dt;

  C1 = 0.65 - 0.78 * rho_local + 0.6 * rho_local * rho_local;
  C2 = 1.06 * rho_local;
  C3 = 1. + 2.6 / sqrt(soil.clay[i]);
  C4 = 0.03 + 0.1 * rho_local * rho_local;
  C5 = 4.;
  soil.k_conductivity_soil[i] = (C1 + 
					C2 * soil.theta[i][1]*(1.-soil.gravel[i]/100.) - 
					(C1 - C4) * exp(-pow(C3*soil.theta[i][1]*(1.-soil.gravel[i]/100.),C5)))
    / (soil.z_soil[i] - soil.z_soil[i-1]);
  /*// ISOLSM formulation
    C1 = (8.80*soil.sand[i]+2.92*soil.clay[i])/(soil.sand[i]+soil.clay[i]);
    C2 = 0.15;
    C3 = (pow(C1,1.-soil.theta_s[i]) * pow(0.6,soil.theta[i][1]) - C2 )
    * soil.theta[i][1]/soil.theta_s[i] + C2;
    soil.k_conductivity_soil[i] = C3  / (soil.z_soil[i] - soil.z_soil[i-1]);*/
    /* printf("LC06 %20.14f\n", rho_local[soil.n_soil]); */

  if (soil.z_litter >= 1e-6)
    {
      // Layer 0 = litter
      //soil.cp_soil[0] = (c_0/rho_0*soil.rho_l + c_w*soil.theta_l[1]) * soil.z_litter / soil.temperature_dt;
      soil.cp_soil[0] = (c_0*soil.rho_l/rho_0 + c_w*soil.theta_l[1]) * soil.z_litter / (2.*soil.temperature_dt);

      // Campbell (1985, Chapter 4)
      C1 = 0.4;
      C2 = 0.5;
      C3 = 1.;
      C4 = 0.06;
      C5 = 4.;
      soil.k_conductivity_soil[0] = (C1 + 
					    C2 * soil.theta_l[1] - 
					    (C1-C4) * exp(-pow(C3*soil.theta_l[1],C5)))
	/ soil.z_litter;
    }

  return;
}


// ----------------------------------------------------
void SOIL_ENERGY_BALANCE()
{
  /*
    The soil energy balance model of Campbell has been adapted to
    compute soil energy fluxes and temperature profiles at the soil
    surface. The model has been converted from BASIC to C. We
    also use an analytical version of the soil surface energy
    balance to solve for LE, H and G.

    Campbell GS (1985) Soil physics with Basic: Transport models for soil-plant systems.
    Elsevier, New York, p 150

    Combine surface energy balance calculations with soil heat
    transfer model to calculate soil conductive and convective heat
    transfer and evaporation rates. Here, only the deep temperature
    is needed and G, Hs and LEs can be derived from air temperature
    and energy inputs.
    Soil evaporation models by Kondo, Mafouf et al. and
    Dammond and Simmonds are used. Dammond and Simmonds for example
    have a convective adjustment to the resistance to heat transfer.
    Our research in Oregon and Canada have shown that this consideration
    is extremely important to compute G and Rn_soil correctly.
  */
  int j, mm1, i, ip1, im1;
  double Fst, Gst;
  double soil_par, soil_nir;
  double u_soil, Rh_soil, kv_soil, kcsoil;
  double a_soil[soilsze+1], b_soil[soilsze+1], c_soil[soilsze+1], d_soil[soilsze+1];
  double est, dest, d2est, tk2, tk3, tk4, llout, lecoef;
  double product;
  double att, btt, ctt;
  double storage;
  double facstab, stabdel, ea;
  double T_new_soil[soilsze+1];
  double k_soil[soilsze+1], cp_soil[soilsze+1], T_soil[soilsze+1];
  double kv_litter, temp;
  double firstlayerevap, secondlayerevap;
  int n_soil_end, k;
//   double temps, templ;

  // Compute soilevap as a function of energy balance at the soil
  // surface. Net incoming short and longwave energy

  // radiation balance at soil in PAR band, W m-2
  //soil_par = (solar.beam_flux_par[1] + solar.par_down[1] - solar.par_up[1]) / 4.6;
  soil_par = (solar.beam_flux_par[1] + solar.par_down[1]) / 4.6
    * (1.-solar.par_soil_refl);

  // radiation balance at soil in NIR band, W m-2
  //soil_nir = (solar.beam_flux_nir[1] + solar.nir_dn[1] - solar.nir_up[1]);
  soil_nir = (solar.beam_flux_nir[1] + solar.nir_dn[1])
    * (1.-solar.nir_soil_refl);

  // incoming radiation balance at soil, solar and terrestrial, W m-2
  soil.rnet = soil_par + soil_nir + solar.ir_dn[1]*epsoil;

  // set air temperature over soil with lowest air layer, filtered
  soil.T_air = prof.tair_filter[1];

  // Compute Rh_soil and rv_soil from wind log profile for lowest layer
  prof.u[1] = UZ(prof.ht[1]/2.);
  u_soil = prof.u[1]; // wind speed one layer above soil

  // Stability factor from Daamen and Simmonds (1996)
  stabdel = 5.*Gravity*(prof.ht[1]/2.)*(soil.tsfc_filter-soil.T_air)/((soil.T_air+TN0)*u_soil*u_soil);
  if (stabdel > 0)
    facstab = 1./pow(1.+stabdel, 0.75);
  else
    facstab = 1./(1.+stabdel)/(1.+stabdel);
  //  facstab = 1./pow(1.+stabdel, 2);
  if (time_var.count <= 1) facstab = 1.;
  facstab = __min(__max(facstab, 0.1), 5.);

  Rh_soil = pow(log((prof.ht[1]/2.-soil.d)/soil.z0), 2.) / vonKarman/vonKarman / u_soil * facstab;
  Rh_soil = __min(__max(Rh_soil, 50.), 3000.);

  // kcsoil is the convective transfer coeff for the soil. (W m-2 K-1)
  kcsoil = (cp * met.air_density) / Rh_soil;

  // soil boundary layer resistance
  soil.rb = Rh_soil;

  // soil resistance
  soil.rs = SOIL_RESISTANCE();

  // litter resistance
  soil.rl = LITTER_RESISTANCE();

  // kvsoil is soil surface conductance to water vapor transfer (s m-1)
  kv_soil = 1. / (soil.rb+soil.rs);

  // kvsoil is litter surface conductance to water vapor transfer (s m-1)
  kv_litter = 1. / (soil.rb+soil.rl);

  // soil relative humidty
  soil.rh_soil = SOIL_RH();

  // litter relative humidty
  soil.rh_litter = LITTER_RH();

  // initialise soil temperatures
  // Note: the first internal soil layer becomes the litter layer
  // all subsequent soil layers are shifted one layer down
  if (soil.z_litter < 1e-6)
    {
      n_soil_end = soil.n_soil;
      for (i=1; i<=soil.n_soil; i++) // does not include litter
        {
          T_soil[i] = soil.T_soil_filter[i];
          k_soil[i] = soil.k_conductivity_soil[i];
          cp_soil[i] = soil.cp_soil[i];
        }
    }
  else
    {
      n_soil_end = soil.n_soil+1;
      for (i=0; i<=soil.n_soil; i++) // includes litter [0] at soil.T_soil, [1] at T_soil
        {
          T_soil[i+1] = soil.T_soil_filter[i];
          k_soil[i+1] = soil.k_conductivity_soil[i];
          cp_soil[i+1] = soil.cp_soil[i];
        }
      T_soil[1] = soil.T_l_filter-TN0; // redundant
    }
  T_soil[0] = soil.T_air;
  T_new_soil[0] = soil.T_air;
  T_soil[n_soil_end+1] = soil.T_base;
  T_new_soil[n_soil_end+1] = soil.T_base;

  // Boundary layer conductance at soil surface, W m-2 K-1
  k_soil[0] = kcsoil;

  // initialize absolute soil surface temperature
  soil.T_Kelvin = soil.T_air+TN0;

  // evaluate latent heat of evaporation at T_soil=T_air
  soil.latent = LAMBDA(soil.T_Kelvin);
  //soil.latent = LAMBDA(soil.tsfc_filter+TN0);

  // evaluate saturation vapor pressure at T_air
  est = ES(soil.T_Kelvin) * 100.; // 100 converts es(T) from mb to Pa

  // Slope of the vapor pressure-temperature curve, Pa/C, evaluate as function of Tk=T_air
  fact.latent = soil.latent;
  dest = DESDT(soil.T_Kelvin);

  // Second derivative of the vapor pressure-temperature curve, Pa/C, evaluate as function of Tk=T_air
  d2est = DES2DT(soil.T_Kelvin);

  // Atmospheric vapour pressure [Pa]
  ea = prof.rhov_air_filter[1][1] * soil.T_Kelvin * Rw;

  // output
  output.c8 = est - ea; // Vapor pressure deficit, Pa

  // Compute products of absolute air temperature
  tk2 = soil.T_Kelvin * soil.T_Kelvin;
  tk3 = tk2 * soil.T_Kelvin;
  tk4 = tk3 * soil.T_Kelvin;

  // Longwave emission at air temperature, W m-2
  llout = epsoil * sigma * tk4;

  // IR emissive flux density from soil [W m-2]
  soil.lout = llout;

  // coefficients for latent heat flux density
  lecoef = met.air_density * 0.622 * soil.latent / met.press_Pa;

  // Weighting factors for solving diff eq.
  Fst = 0.6;
  Gst = 1. - Fst;

  // solve by looping through the d[]/dt term of the Fourier
  // heat transfer equation

  soil.soilevap   = soil.soilevap_filter;
  soil.litterevap = soil.litterevap_filter;
  for (j=1; j<=soil.temperature_mtime; j++)
    {
      if (soil.z_litter < 1e-6)
        {
          firstlayerevap  = soil.soilevap;
          secondlayerevap = 0.;
        }
      else
        {
          firstlayerevap  = soil.litterevap;
          secondlayerevap = soil.soilevap;
        }
      for (i=1; i<=n_soil_end; i++)
        {
          // define coef for each soil layer
          im1 = i-1;
          ip1 = i+1;
          c_soil[i]   = -k_soil[i]*Fst;
          a_soil[ip1] = c_soil[i];
          b_soil[i]   = Fst*(k_soil[i]+k_soil[im1]) + cp_soil[i];
          d_soil[i]   = Gst*k_soil[im1]*T_soil[im1] + (cp_soil[i]-Gst*(k_soil[i] + k_soil[im1]))*T_soil[i]
            + Gst*k_soil[i]*T_soil[ip1];
        }
      d_soil[1]          = d_soil[1] + k_soil[0]*T_new_soil[0]*Fst + soil.rnet - soil.lout - firstlayerevap;
      d_soil[2]          = d_soil[2] - secondlayerevap;
      d_soil[n_soil_end] = d_soil[n_soil_end] + k_soil[n_soil_end]*Fst*T_new_soil[n_soil_end+1];

      mm1 = n_soil_end-1;
      for (i=1; i<=mm1; i++)
        {
          ip1 = i+1;
          c_soil[i]   = c_soil[i]/b_soil[i];
          d_soil[i]   = d_soil[i]/b_soil[i];
          b_soil[ip1] = b_soil[ip1] - a_soil[ip1]*c_soil[i];
          d_soil[ip1] = d_soil[ip1] - a_soil[ip1]*d_soil[i];
        }
      T_new_soil[n_soil_end] = d_soil[n_soil_end]/b_soil[n_soil_end];

      for (i=mm1; i>=1; i--)
        T_new_soil[i] = d_soil[i] - c_soil[i]*T_new_soil[i+1];

      // soil temperature at 15 cm
      soil.T_15cm = T_new_soil[7];

      // compute soil conductive heat flux density, W m-2
      if (soil.z_litter < 1e-6)
	{
	  soil.gsoil = k_soil[1]*(T_new_soil[1]-T_new_soil[2]); // Matthias ??? why 1 and 2
	  // Does not work without litter, do not know why
	  //storage = 0.;
	  //for (k=1; k<=n_soil_end; k++)
	  // storage += cp_soil[k]*(T_new_soil[k]-T_soil[k]);
	  //soil.gsoil += storage;
	}
      else
	{ // take 2 and 3 because 1 is litter
	  soil.gsoil = k_soil[2]*(T_new_soil[2]-T_new_soil[3]); 
	  storage = 0.;
	  for (k=2; k<=n_soil_end; k++)
	    storage += cp_soil[k]*(T_new_soil[k]-T_soil[k]);
	  soil.gsoil += storage;
	}

      //output soil heat flux [W m-1]
      output.c7 = soil.gsoil;

      // solve for T_litter using quadratic solution x1/2 = (-b +- sqrt(b^2-4ac) )/2a

      att = 6.*epsoil*sigma*tk2 + d2est/2.*lecoef*(kv_soil*soil.rh_soil+soil.c_litterevap*kv_litter*soil.rh_litter);

      btt = 4.*epsoil*sigma*tk3 + kcsoil +
        lecoef*(dest*(kv_soil*soil.rh_soil+soil.c_litterevap*kv_litter*soil.rh_litter)-kv_soil*soil.rh_soil*d2est*soil.gsoil/k_soil[1]);

//       ctt = llout - soil.rnet + soil.gsoil
//         + lecoef*(__max(soil.rh_soil*est-ea, 0.)*kv_soil+soil.c_litterevap*__max(soil.rh_litter*est-ea, 0.)*kv_litter)
//         + lecoef*kv_soil*soil.rh_soil*(d2est/2.*pow(soil.gsoil/k_soil[1], 2) - dest*soil.gsoil/k_soil[1]);
      ctt = llout - soil.rnet + soil.gsoil
        + lecoef*((soil.rh_soil*est-ea)*kv_soil+soil.c_litterevap*(soil.rh_litter*est-ea)*kv_litter)
        + lecoef*kv_soil*soil.rh_soil*(d2est/2.*soil.gsoil/k_soil[1]*soil.gsoil/k_soil[1]
					      - dest*soil.gsoil/k_soil[1]);

      product = btt*btt - 4.*att*ctt;

      if (product >= 0.)
        soil.tsfc = soil.T_air + (-btt+sqrt(product))/(2.*att);
      else
        soil.tsfc = soil.T_air;

      /* ToDo for isotopes // litter evaporation from second Taylor expansion
      soil.litterevap = soil.c_litterevap*lecoef*kv_litter*(__max(soil.rh_litter*est-ea, 0.)
                                                                          + soil.rh_litter*dest*(soil.tsfc-soil.T_air)
                                                                          + soil.rh_litter*d2est/2.*pow(soil.tsfc-soil.T_air, 2)
                                                                          ); */
      // litter evaporation from energy budget
      if (1 == 1)
	soil.litterevap = soil.c_litterevap*lecoef*kv_litter
	  * (soil.rh_litter*ES(soil.tsfc+TN0)*100.
	     - prof.rhov_air_filter[1][1]*(soil.T_air+TN0)*Rw);
      else
	soil.litterevap = soil.c_litterevap*lecoef*kv_litter
	  * ((soil.rh_litter*est-ea)
	     + soil.rh_litter*dest*(soil.tsfc-soil.T_air)
	     + soil.rh_litter*d2est/2.*pow(soil.tsfc-soil.T_air, 2));

      if (set_switch.no_neg_water_flux == 1)
	soil.litterevap = __max(soil.litterevap, 0.);
  
      /* ToDo for isotopes // soil evaporation from second Taylor expansion
      soil.soilevap = lecoef*kv_soil*(__max(soil.rh_soil*est-ea, 0.)
                                             + soil.rh_soil*dest*(temp-soil.T_air)
                                             + soil.rh_soil*d2est/2.*pow(temp-soil.T_air, 2)
                                             ); */
      // soil evaporation from energy budget
      if (1 == 1)
	soil.soilevap = lecoef*kv_soil
	  * (soil.rh_soil*ES(soil.tsfc+TN0)*100.
	     - prof.rhov_air_filter[1][1]*(soil.T_air+TN0)*Rw);
      else
	{
	  temp = soil.tsfc - soil.gsoil/k_soil[1];
	  soil.soilevap = lecoef*kv_soil
	    * ((soil.rh_soil*est-ea)
	       + soil.rh_soil*dest*(temp-soil.T_air)
	       + soil.rh_soil*d2est/2.*pow(temp-soil.T_air, 2));
	}
      if (set_switch.no_neg_water_flux == 1)
	soil.soilevap = __max(soil.soilevap, 0.);

      soil.evap = soil.litterevap + soil.soilevap;

      // IR emissive flux density from soil [W m-2]
      //soil.lout = epsoil * sigma * pow(soil.tsfc+TN0, 4);
      soil.lout = llout
        + 4.*epsoil*sigma*tk3*(soil.tsfc-soil.T_air)
        + 6.*epsoil*sigma*tk2
	*(soil.tsfc-soil.T_air)*(soil.tsfc-soil.T_air);

      // Sensible heat flux density over soil [W m-2]
      soil.heat = kcsoil*(soil.tsfc-soil.T_air);
      
      // update latent heat
      //soil.latent = LAMBDA(soil.tsfc_filter+TN0);
      //lecoef = met.air_density * 0.622 * soil.latent / met.press_Pa;

      // set new temperature profile, check for extremes
      for (i=0; i<=n_soil_end+1; i++)
        {
          if (T_new_soil[i] < -30. || T_new_soil[i] > 60.)
	    T_new_soil[i] = input.ta;
          T_soil[i] = T_new_soil[i];
        }
    } // next j of soil.temperature_mtime
  
  // set 'real' soil & litter temperatures
  if (soil.z_litter < 1e-6)
    {
      soil.T_l = soil.tsfc+TN0;
      for (i=0; i<=n_soil_end; i++)
	soil.T_soil[i] = T_new_soil[i+1];
    }
  else
    {
      soil.T_l = T_new_soil[1]+TN0;
      for (i=0; i<=n_soil_end; i++)
	soil.T_soil[i] = T_new_soil[i+1];
    }

  // store coefficient of latent heat
  soil.lecoef = lecoef;

  // for output
  output.c1 = kcsoil;
  output.c2 = kv_soil;
  output.c3 = Rh_soil;
  output.c4 = soil.rb;
  output.c5 = met.air_density;

  // Restrict litterevap to maximum available litter water content
  if (soil.z_litter < 1e-6)
    soil.c_litterevap = soil.c_litterevap_save = 0.;
  else
    {
      soil.c_litterevap_save = soil.c_litterevap;
      if (soil.litterevap > 0.)
	{
	  //* soil.latent; // kg/m2/s = mm/s -> W/m2
	  soil.maxlitterevap = (1000.*soil.z_litter/time_var.time_step * soil.theta_l[1] // water in litter
				       - __min(__max(0.,                                                   // - drainage in soil
						     1000.*soil.z_litter/time_var.time_step*
						     (soil.theta_l[1]-soil.theta_l33)),
					       soil.qinfl_max))
	    * soil.latent; // kg/m2/s = mm/s -> W/m2
	  if (soil.c_litterevap != 0.)
	    temp = soil.litterevap / soil.c_litterevap;
	  else
	    temp = 0.;
	  soil.c_litterevap = __min(1., soil.maxlitterevap/__max(temp, 1e-20));
	  soil.maxlitter = 0;
	  if (soil.litterevap > soil.maxlitterevap)
	    {
	      soil.litterevap_save = soil.litterevap;
	      soil.litterevap = soil.maxlitterevap;
	      soil.maxlitter = 1;
#ifdef DIAG
	      printf("Litterevap %e > Maxlitterevap %e at %ld.", soil.litterevap, 
		     soil.maxlitterevap, time_var.daytime);
#endif
	    }
	}
    }

  // W/m2 -> kg/m2/s = mm/s
  flux.soilevap[1] = soil.soilevap / soil.latent;
  flux.litterevap[1] = soil.litterevap / soil.latent;
  flux.s_evap[1] = flux.litterevap[1] + flux.soilevap[1];

  return;
}


// ----------------------------------------------------
double UZ(double zzz)
{
  // Wind speed inside canopy is calculated based on Campbell and Norman, 1998
  double y, zh, factor, uh, aa, ulow, ustarlow;

  zh = zzz/ht;

  if (zh >= 0.11)
    {
      factor = log((ht-zd)/z0)/vonKarman;

      uh = factor*met.ustar; // wind speed at canopy top

      // calculating attenuation factor results in very large values (>8) due to the height canopy.
      // Measured values are typically around 2.5

      aa = (attfac-1.5)*time_var.lai/lai+1.5; // vary attenuation factor between 1.5 and attfac (from parameterfile)

      y = uh * exp(aa*(zh-1)); // calculate wind speed at height zh with attenuation
    }
  else
    {
      ulow = UZ(0.11*ht); // calculate wind speed at 0.11 * canopy height
      ustarlow = ulow*vonKarman/log((0.11*ht/0.01)); // calculate ustar at 0.11*canopy height
      y = ustarlow/vonKarman*log((zh*ht)/0.01); // calculate wind speed at height zh with logarithmic wind profile
    }

  // keep routines from blowing up
  if (y < 0.1)
    y = 0.1;

  return y;
}


// ----------------------------------------------------
void CARBON_ISOTOPE()
{
  /*
  This subroutine, ISOTOPE, computes concentrations of 13C
  It was developed in collaboration with David Bowling (circa Dec, 2000).

  This program computes fluxes as 13C/(13C+12C)

  We apply the algorithm of Farquhar et al that computes D13 as a function of Ci/Ca.
  Ci/Ca is computed using the photosynthesis, respiration and leaf energy balance algorithms
  of CANOAK. The Farquhar photosynthesis model is coupled to the Ball-Berry stomatal conductance
  model.

  A Lagrangian random walk model is used to compute the turbulence Dispersion matrix
  and compute concentration profiles in and above the canopy. This scheme accounts
  for counter gradient transfer.
   
  Summary of variables computed:
    Isotopic discrimination due to photosynthesis and diffusion
    D13(big Delta) = a + (b-a) ci/ca
    a is fractionation by diffusion, 4.4 per mill
    b is net fractionation by carboxylation, typically 27 per mill
    Isotopic discrimination in term of the isotopic content of the air and plant
    D13 = (d13_air - d13_plant)/(1+d13_plant)
    D13 = (Rair/Rplant-1) 1000
    D13 = alpha -1
    Rplant = Rair/(D13+1)
    Isotopic content relative to the PeeDee Standard
    d13_plant(little delta)= (Rplant/Rpdb-1) * 1000 (per mill)
    d13_air(little delta)= (Rair/Rpdb-1) * 1000 (per mill)
    Rpdb_13 = 0.01124 (13C/12C; Farquhar et al. 1989) or
    Rpdb_12_13 = 0.01115 (13C/(12C+13C), Tans et al)
    The ratio Rpdb_13/Rpdb_12_13 = 1.00807. It can be used as a multiplier to convert the
      isotopic ratio.
    alpha= Rair/Rplant
    photosynthetic flux of 13C= A*Ra/(1+D13) = A * Rplant
      This relation is in terms of net Ps, A, not gross Ps
  */

  double R_soil_12C, R_air_12C, co2_13; //, da13, R_soil_CO2, R_air_CO2;
  double b_a, al, as, ab, a, b, source_sun=0., source_shade=0.;
  double sun_A=0, shd_A=0., sumA=0;
  int j, timelag;

  /*
    R_soil_C = 13C/(12C) is the isotopic ratio of the soil C flux and is calculated as
    the sum of heterotropic signal (constant) and autotrophic signal from prevoius day with time lag (=canopy discrimination)
  */
  timelag = 1;  // for now we set timelag to 1 day,but depending on site could be longer, Hainich=4 days

  R_soil_12C = INVDELTA1000(Cisotope.delta_soil)*Rpdb_12C*soil.respiration_hetero +
    INVDELTA1000(Cisotope.bigdelta_long[timelag])*Rpdb_12C*soil.respiration_auto;

  // soil respiration is for total CO2, so we need to use the
  // Tans version of Pee Dee Bee, 13C/(12C + 13C) to compute the flux
  // of 13 C for respiration

  //R_soil_CO2 = (Cisotope.delta_soil/1000.+1.)*Rpdb_CO2;

  //R_soil_CO2 = R_soil_C12/(1+R_soil_C12);


  /*
    co2_13 is upper boundary condition for [13CO2], in ppm
  */

  R_air_12C = INVDELTA1000(input.d13CO2)*Rpdb_12C; // ratio of 13C relative to 12C
  //R_air_CO2 = R_air_C12/(1+R_air_C12);  // convert to 13C relatve to (12C + 13C)

  /*
    Or calculate it with Bowling method
    Calculate calculate d13Cair from regression of deltaCa*Ca=m*Ca+b
    and then solve for delta
    Data are coming from field measurements and coefficients are provided in parameter file
    Note different functions exist for day and night

  if(input.parin > 0)
      da13 = (Cisotope.da_m_day * input.co2air + Cisotope.da_b_day)/input.co2air;  // daytime
  else
      da13= (Cisotope.da_m_night * input.co2air + Cisotope.da_b_night)/input.co2air;   //  nighttime

  R_air_12C = (da13/1000.+1.)*Rpdb_12C;  // ratio of 13C relative to CO2
  */


  /*
    We have to be careful here as CO2air is 12C + 13C. We need to multiply
    R_13_12 times 12C to compute the concentration of [13CO2]

    R_13_12 * 12C = 13C
    C = 12C + 13C
    12C = C - 13C
    12C= C/(1+R_13_12)
    13C = R_13_12 * 12C = R_13_12 * (C - 13C)
    13C + R_13_12 * 13C = R_13_12 * C
    13C = R_13_12 * C/(1+ R_13_12)
  */


  co2_13 = R_air_12C * input.co2air; // ppm concentration of 13C

  // Respiration flux of 13C = Rr*Resp, where Resp = soil respiration

  soil.resp_13 = R_soil_12C * soil.respiration_mole;

  
  b_a=(27-4.4);
  ab = 2.9; //fractionation due to diffusion in the laminar boundary layer, permille
  a = 4.4; //fractionation due to diffusion from the leaves to stomata cavity, permille
  as = 1.1; //fractionation as CO2 enters solution at 25 deg C
  al = 0.7; //fractionation by diffusion within cell
  b = 29.; //fractionation during carboxylation

  // calculate canopy discrimination (only if there is photosynthesis)

  //  if (flux.photosyn > 0.) {

  // compute cc (CO2 mixing ratio at site of carboxylation) on sun and shade fractions for layers of the canopy
  // and use this information to compute the source/sink function of 13C

  for (j=1; j<=jtot; j++)
    {
      /*
      // get resistances and convert from m s-1 to mol m-2 s-1
      fact_rb = rugc * (prof.tair[j]+TN0)/ met.press_Pa; //conversion factor
      // J mol-1 K-1 * K * Pa-1 = J mol-1 Pa-1 = (Kg m2 s-2) mol-1 (Kg m-1 s-2)-1 = m3 mol-1
      fact_rs_sun = rugc * (prof.sun_tleaf_filter[j]+TN0)/ met.press_Pa; //conversion factor
      fact_rs_shd = rugc * (prof.shd_tleaf_filter[j]+TN0)/ met.press_Pa; //conversion factor
      rb = prof.sun_rbco2[j]*fact_rs_sun; // s m-1 to mol-1 m2 s1
      rs = prof.sun_rs[j]*fact_rs_sun*1.577; // s m-1 to mol-1 m2 s1 and convert from rs(H2O) to rs(CO2)
      cs = prof.co2_air[j] - rb*prof.sun_A[j]; //unit ppm
      prof.sun_ci[j] = prof.shd_cs[j] - rs*prof.sun_A[j]; //unit ppm
      ci = cs - rs*prof.sun_A[j]; //unit ppm
      */

      /* DONE IN PHOTOSYNTHESIS NOW
      // mesophyll conductance
      rm = 1/(0.0045*prof.vcmax[j]); // mol-1 m2 s1 similar to Suits et al. 2005
      rm_epron = 80/(prof.sun_A[j]-0.05); // according to Epron et al. 1995
      rm_manter = 1/(0.0028*prof.vcmax[j] - 0.0176); // accornding to data from Manter 2004;

      fact_rs_sun = rugc * (prof.sun_tleaf_filter[j]+TN0)/ met.press_Pa; //conversion factor
      fact_rs_shd = rugc * (prof.shd_tleaf_filter[j]+TN0)/ met.press_Pa; //conversion factor

      rm_gs_sun = (1/0.9)*prof.sun_rs[j]*fact_rs_sun*1.577; // s m-1 to mol-1 m2 s1 and convert from rs(H2O) to rs(CO2), according to Manter (data, pers. comm)
      rm_gs_shd = (1/0.9)*prof.shd_rs[j]*fact_rs_shd*1.577; // s m-1 to mol-1 m2 s1 and convert from rs(H2O) to rs(CO2), according to Manter (data, pers. comm)

      rm_sun = rm;
      rm_shd = rm;

      //sun leaves
      prof.sun_cc[j] = prof.sun_ci[j] - rm_sun*prof.sun_A[j]; //unit ppm

      //shade leaves
      prof.shd_cc[j] = prof.shd_ci[j] - rm_shd*prof.shd_A[j]; //unit ppm

      // Cc/Ca on sun leaves
      prof.sun_cica[j]= prof.sun_ci[j]/prof.co2_air[j];
      */

      // Compute Del 13C using Farquhar eq. D13C= a + (b-a)Ci/Ca
      // a = 4.4 per mill, b = 27.5 per mill

      prof.sun_D13[j] = a + b_a*__min(prof.sun_cica[j], 1); // parts per mill

      // compute Del 13C using Farquhar eq Appendix

      prof.sun_D13_ab[j] = ab*__max((prof.co2_air_filter[j]-prof.sun_cs[j])/prof.co2_air_filter[j], 0.);
      prof.sun_D13_a[j] = a*__max((prof.sun_cs[j]-prof.sun_ci[j])/prof.co2_air_filter[j], 0.);
      prof.sun_D13_asal[j] = (as+al)*__max((prof.sun_ci[j]-prof.sun_cc[j])/prof.co2_air_filter[j], 0.);
      prof.sun_D13_b[j] = b*__min((prof.sun_cc[j])/prof.co2_air_filter[j], 1.);

      prof.sun_D13_long[j] = prof.sun_D13_ab[j] + prof.sun_D13_a[j] + prof.sun_D13_asal[j] + prof.sun_D13_b[j];

      // Compute Del 13C using Farquhar eq. D13C= a + (b-a)Ci/Ca

      prof.shd_D13[j]=a+b_a*__min(prof.shd_cica[j], 1.);

      prof.shd_D13_ab[j] = ab*__max((prof.co2_air_filter[j]-prof.shd_cs[j])/prof.co2_air_filter[j], 0.);
      prof.shd_D13_a[j] = a*__max((prof.shd_cs[j]-prof.shd_ci[j])/prof.co2_air_filter[j], 0.);
      prof.shd_D13_asal[j] = (as+al)*__max((prof.shd_ci[j]-prof.shd_cc[j])/prof.co2_air_filter[j], 0.);
      prof.shd_D13_b[j] = b*__min((prof.shd_cc[j])/prof.co2_air_filter[j], 1.);

      prof.shd_D13_long[j] = prof.shd_D13_ab[j] + prof.shd_D13_a[j] + prof.shd_D13_asal[j] + prof.shd_D13_b[j];

      /*
        sour13co2 is photsynthetic flux of 13C
        D13 is integrated (sun+shade) discrimination for this layer
        need to work up D13 from prof.sun_D13 and prof.sh_D13!
        - may need to iterate between sour13CO2 and Rair
      */

      // need to divide D13 by 1000

      prof.Rplant_sun[j] =prof.R13_12_air[j]/(1.+ prof.sun_D13_long[j]/1000.);

      prof.Rplant_shd[j] =prof.R13_12_air[j]/(1.+ prof.shd_D13_long[j]/1000.);


      // sour13CO2 source/sink strength of 13C in micromoles
      // apply minus sign as Ps is an atmospheric sink for C

      // fluxes on the sun/shade fractions are weighted and considered

      sun_A = __max(prof.sun_A[j]*solar.prob_beam[j], 0.);
      shd_A = __max(prof.shd_A[j]*solar.prob_shd[j], 0.);
      sumA = sun_A + shd_A;
      if (sun_A > 0. || shd_A > 0.)
        {
          prof.D13C_long[j]=(shd_A*prof.shd_D13_long[j] + sun_A*prof.sun_D13_long[j]) / sumA;
          prof.D13C[j]=(shd_A*prof.shd_D13[j] + sun_A*prof.sun_D13[j]) / sumA;

          // step of discrimination
          prof.D13C_ab[j]=(shd_A*prof.shd_D13_ab[j] + sun_A*prof.sun_D13_ab[j]) / sumA;
          prof.D13C_a[j]=(shd_A*prof.shd_D13_a[j] + sun_A*prof.sun_D13_a[j]) / sumA;
          prof.D13C_asal[j]=(shd_A*prof.shd_D13_asal[j] + sun_A*prof.sun_D13_asal[j]) / sumA;
          prof.D13C_b[j]=(shd_A*prof.shd_D13_b[j] + sun_A*prof.sun_D13_b[j]) / sumA;

          // cs, ci, cc
          prof.cs[j]=(shd_A*prof.shd_cs[j] + sun_A*prof.sun_cs[j]) / sumA;
          prof.ci[j]=(shd_A*prof.shd_ci[j] + sun_A*prof.sun_ci[j]) / sumA;
          prof.cc[j]=(shd_A*prof.shd_cc[j] + sun_A*prof.sun_cc[j]) / sumA;

	  prof.csca[j] = prof.cs[j]/prof.co2_air_filter[j];
	  prof.cica[j] = prof.ci[j]/prof.co2_air_filter[j];
	  prof.ccca[j] = prof.cc[j]/prof.co2_air_filter[j];

          source_shade= shd_A*prof.Rplant_shd[j];
          source_sun = sun_A*prof.Rplant_sun[j];
        }
      else
	{
	  prof.D13C_long[j] = undef;
	  prof.D13C[j] = undef;
	  // step of discrimination
	  prof.D13C_ab[j] = undef;
	  prof.D13C_a[j] = undef;
	  prof.D13C_asal[j] = undef;
	  prof.D13C_b[j] = undef;
	  // cs, ci, cc
	  prof.cs[j] = undef;
	  prof.ci[j] = undef;
	  prof.cc[j] = undef;
	  prof.csca[j] = undef;
	  prof.cica[j] = undef; 
	  prof.ccca[j] = undef; 
	}

      // if dark, respire only fixed 13C; assume signature of soil for now
      // latter use data assimilated from day

      if (prof.shd_A[j] <= 0.)
        source_shade = solar.prob_shd[j] * prof.shd_A[j] * prof.Rresp_ave[j];

      if (prof.sun_A[j] <= 0.)
        source_sun = solar.prob_beam[j] * prof.sun_A[j] * prof.Rresp_ave[j];

      // micromol C m-2 s-1

      prof.sour13co2[j] = -prof.dLAIdz[j]*(source_shade+source_sun);
    }

  //  } // end if photosyntheis


  // add bole respiration to 13CO2 sources assuming it has the canopy discrimination signal
  // shifted by one day
  
  for (j=1; j<=jtot; j++)
    prof.sour13co2[j] += bole.layer[j]*INVDELTA1000(Cisotope.bigdelta_long[timelag])*Rpdb_12C;

  // generate 13CO2 concentration profiles using information on 13C leaf
  // sources and sinks, turbulent mixing of the air and source of soil
  // respiration, each with distinct signatures


  // need to convert umol m-3 to umol/mol we have to consider
  // Pc/Pa = [CO2]ppm = rhoc ma/ rhoa mc

  fact.co2 = (mass_air/mass_13CO2)*met.air_density_mole;
  CONC(prof.sour13co2, prof.c13cnc, co2_13, soil.resp_13, fact.co2);

  for (j=1; j<=jtot3; j++)
    {
      // this is the 13C/(12C+13C) ratio of air in the canopy
      prof.R13_12_air[j] = prof.c13cnc[j]/prof.co2_air_filter[j];
      prof.d13Cair[j] = (prof.R13_12_air[j]/Rpdb_CO2-1.)*1000.;
    }
  /* printf("d13a: %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f\n", prof.R13_12_air[1], prof.c13cnc[1], prof.co2_air_filter[1], prof.sour13co2[1], bole.layer[1], INVDELTA1000(Cisotope.bigdelta_long[timelag])*Rpdb_12C); */

  return;
}


// ----------------------------------------------------
void READ_PARAMETER()
{
  int param_counter;
  char record[256], *field;
  FILE *param_file;

  param_file = fopen(parameter_file,"r");
  if (param_file == NULL) puts("Cannot open parameter file.txt.\n");

  param_counter = 1;
  while (!feof(param_file))
    {
      while (fgets(record, sizeof(record), param_file))
        {
          if (ferror(param_file))
            {
              perror("Error reading parameter file");
              exit(1);
              break;
            }
          if (record[0] != '#' && record[0] != ';' &&
	      record[0] != ',' && record[0] != '\0' &&
	      record[0] != '\n')
            {
              if (record[strlen(record) - 1] != '\n') /* make sure the newline is there */
                record[strlen(record) - 1] = '\n';
              field = strtok(record, "=");
              if (strlen(field) != 0)
                {
		  strncpy(parameter_name[param_counter], 
			  ltrim(field),
			  (size_t) find_last_nonwhite_character(field)+1);
		  field = strtok(NULL, ";");
		  if (field == NULL)
		    {
		      printf("Problem reading parameter.\n");
		      printf("  Counter: %i\n", param_counter);
		      printf("  Name: %s\n", parameter_name[param_counter]);
		      printf("  Line: %s\n", record);
		    }
		  else
		    {
		      // correct last character if no ; found
		      if (field[strlen(field)-1] == '\n')
			field[strlen(field)-1] = '\0';
		      strncpy(parameter_value[param_counter], 
			      ltrim(field),
			      (size_t) find_last_nonwhite_character(field)+1);
		      param_counter++;
		    }
                }
            }
        }
    }

  param_counter--;
  nparameter = param_counter;

  return;
}


// ----------------------------------------------------
void SET_PARAMETER()
{
  // ------------------------------------------------------
  // DECLARE and read PARAMETER VALUES
  // ------------------------------------------------------
  //
  // Parameter that were read in from parameter_file in READ_PARAMETER()
  // are assigned to internal parameters overwriting the default values.
  int i, j;
  char stmp[250];
  double vc25_up=0., vc25_down=0., jm_vc_up=0., jm_vc_down=0., rd_vc_up=0., rd_vc_down=0.;
  double g0_up=0., g0_down=0., a1_up=0., a1_down=0., D0_up=0., D0_down=0.;
  double kball_up=0., kball_down=0., bprime_up=0., bprime_down=0.;

  // first loop through parameters
  for (i=1; i<=nparameter; i++)
    {
      // directories & suffixes
      if (strcmp(parameter_name[i], "workdir") == 0)
	strcpy(disksubdir, parameter_value[i]);
      if (strcmp(parameter_name[i], "outsuffix") == 0)
	strcpy(filesuffix, parameter_value[i]);
      if (strcmp(parameter_name[i], "indir") == 0)
	strcpy(input_file_subdir, parameter_value[i]);
      if (strcmp(parameter_name[i], "insuffix") == 0)
	strcpy(input_file_suffix, parameter_value[i]);
      if (strcmp(parameter_name[i], "dispfile") == 0)
	strcpy(dispfile, parameter_value[i]);
      // model set-up
      if (strcmp(parameter_name[i], "year0") == 0)
	time_var.year0 = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "start_run") == 0)
	start_run = atol(parameter_value[i]);
      if (strcmp(parameter_name[i], "end_run") == 0)
	end_run = atol(parameter_value[i]);
      if (strcmp(parameter_name[i], "start_profiles") == 0)
	start_profiles = atol(parameter_value[i]);
      if (strcmp(parameter_name[i], "end_profiles") == 0)
	end_profiles = atol(parameter_value[i]);
      if (strcmp(parameter_name[i], "time_step") == 0)
	time_var.time_step = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "switch_soil_resp_temp") == 0)
	set_switch.soil_resp_temp = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "switch_mode") == 0)
	set_switch.mode = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "switch_ball") == 0)
	set_switch.ball = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "switch_isoprene") == 0)
	set_switch.isoprene = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "switch_d13c") == 0)
	set_switch.d13c = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "switch_wiso") == 0)
	set_switch.wiso = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "switch_bethy_resp") == 0)
	set_switch.bethy_resp = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "switch_no_negative_water_flux") == 0)
	set_switch.no_neg_water_flux = atoi(parameter_value[i]);
      // ecosystem - basic
      if (strcmp(parameter_name[i], "latitude") == 0)
	latitude = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "longitude") == 0)
	longitude = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "zone") == 0)
	zone = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "ht") == 0)
	ht = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "pai") == 0)
	pai = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "lai") == 0)
	lai = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "vc25") == 0)
	vc25[1] = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "jm_vc") == 0)
	jm_vc[1] = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "rd_vc") == 0)
	rd_vc[1] = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "ustar_ref") == 0)
	ustar_ref = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "jtot") == 0)
	jtot = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "zm") == 0)
	zm = (double) atof(parameter_value[i]);
      // ecosystem - photosynthesis
      if (strcmp(parameter_name[i], "hkin") == 0)
	hkin = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "skin") == 0)
	skin = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "ejm") == 0)
	ejm = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "evc") == 0)
	evc = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "kc25") == 0)
	kc25 = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "ko25") == 0)
	ko25 = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "o2") == 0)
	o2 = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "tau25") == 0)
	tau25 = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "ekc") == 0)
	ekc = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "eko") == 0)
	eko = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "erd") == 0)
	erd = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "ektau") == 0)
	ektau = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "tk_25") == 0)
	tk_25 = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "toptvc") == 0)
	toptvc = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "toptjm") == 0)
	toptjm = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "curvature") == 0)
	curvature = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "qalpha") == 0)
	qalpha = (double) atof(parameter_value[i]);
      // ecosystem - stomata
      if (strcmp(parameter_name[i], "g0") == 0)
	g0[1] = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "a1") == 0)
	a1[1] = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "D0") == 0)
	D0[1] = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "gm_vc") == 0)
	gm_vc = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "fthreshold") == 0)
	sfc_res.fthreshold = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "fslope") == 0)
	sfc_res.fslope = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "kball") == 0)
	kball[1] = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "bprime") == 0)
	bprime[1] = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "rsm") == 0)
	rsm = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "brs") == 0)
	brs = (double) atof(parameter_value[i]);
      // ecosystem - leaf
      if (strcmp(parameter_name[i], "ep") == 0)
	ep = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "n_stomata_sides") == 0)
	n_stomata_sides = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "betfact") == 0)
	betfact = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "markov") == 0)
	markov = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "lleaf") == 0)
	lleaf = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "leaf_out") == 0)
	leaf_out = atol(parameter_value[i]);
      if (strcmp(parameter_name[i], "leaf_full") == 0)
	leaf_full = atol(parameter_value[i]);
      if (strcmp(parameter_name[i], "leaf_fall") == 0)
	leaf_fall = atol(parameter_value[i]);
      if (strcmp(parameter_name[i], "leaf_fall_complete") == 0)
	leaf_fall_complete = atol(parameter_value[i]);
      // ecosystem - misc
      if (strcmp(parameter_name[i], "attfac") == 0)
	attfac = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "eabole") == 0)
	eabole = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "epsoil") == 0)
	epsoil = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "water_film_thickness") == 0)
	water_film_thickness = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "tau_water") == 0)
	tau_water = (double) atof(parameter_value[i]);
      // 13co2
      if (strcmp(parameter_name[i], "delta_soil") == 0)
	Cisotope.delta_soil = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "da_m_day") == 0)
	Cisotope.da_m_day = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "da_b_day") == 0)
	Cisotope.da_b_day = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "da_m_night") == 0)
	Cisotope.da_m_night = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "da_b_night") == 0)
	Cisotope.da_b_night = (double) atof(parameter_value[i]);
      // litter
      if (strcmp(parameter_name[i], "z_litter") == 0)
	soil.z_litter = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "rho_l") == 0)
	soil.rho_l = (double) atof(parameter_value[i]);
      // soil
      if (strcmp(parameter_name[i], "saxton") == 0)
	soil.saxton = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "n_soil") == 0)
	soil.n_soil = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "root_factor") == 0)
	soil.root_factor = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "clay_factor") == 0)
	soil.clay_factor = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "sand_factor") == 0)
	soil.sand_factor = (double) atof(parameter_value[i]);
      // water isotopes
      if (strcmp(parameter_name[i], "wiso_nofracsoil") == 0)
	wiso.nofracsoil = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "wiso_nofraclitter") == 0)
	wiso.nofraclitter = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "wiso_nofracleaf") == 0)
	wiso.nofracleaf = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "wiso_nofracin") == 0)
	wiso.nofracin = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "wiso_implicit") == 0)
	wiso.implicit = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "merlivat") == 0)
	wiso.merlivat = atoi(parameter_value[i]);
      // Nate McDowell extra
      if (strcmp(parameter_name[i], "extra_nate") == 0)
	extra_nate = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "nup") == 0)
	nup = atoi(parameter_value[i]);
      if (strcmp(parameter_name[i], "vc25_up") == 0)
	vc25_up = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "vc25_down") == 0)
	vc25_down = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "jm_vc_up") == 0)
	jm_vc_up = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "jm_vc_down") == 0)
	jm_vc_down = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "rd_vc_up") == 0)
	rd_vc_up = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "rd_vc_down") == 0)
	rd_vc_down = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "g0_up") == 0)
	g0_up = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "g0_down") == 0)
	g0_down = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "a1_up") == 0)
	a1_up = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "a1_down") == 0)
	a1_down = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "D0_up") == 0)
	D0_up = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "D0_down") == 0)
	D0_down = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "kball_up") == 0)
	kball_up = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "kball_down") == 0)
	kball_down = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "bprime_up") == 0)
	bprime_up = (double) atof(parameter_value[i]);
      if (strcmp(parameter_name[i], "bprime_down") == 0)
	bprime_down = (double) atof(parameter_value[i]);
    }

  // derived variables from first pass through parameters

  zd = 2.0/3.0*ht; // displacement height [m]
  z0 = 0.1*ht; // rougness lenght [m]

  jtot1 = jtot + 1;
  jtot3 = jtot*3; // number of layers in the domain, three times canopy height

  izref = (int) ceil(zm*jtot/ht); // array value of reference height = measurement height*jtot/ht

  //MC20190702 pi4  = 12.5663706;
  delz = ht/jtot; // height of each layer, ht/jtot
  zh65 = 0.65/ht; // 0.65/ht

  epsigma   = ep*sigma; //5.5566e-8; // ep*sigma
  epsigma2  = 2*ep*sigma; //11.1132e-8; // 2*ep*sigma
  epsigma4  = 4*ep*sigma; //22.2264e-8; // 4.0 * ep * sigma
  epsigma6  = 6*ep*sigma; //33.3396e-8; // 6.0 * ep * sigma
  epsigma8  = 8*ep*sigma; //44.448e-8; // 8.0 * ep * sigma
  epsigma12 = 12*ep*sigma; //66.6792e-8; // 12.0 * ep * sigma

  qalpha2 = qalpha*qalpha;

  // second loop through parameters, n_soil times
  soil.z_soil[0] = 0.;
  for (j=1; j<=soil.n_soil; j++)
    {
      for (i=1; i<=nparameter; i++)
	{
	  sprintf(stmp, "z_soil%02i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    soil.z_soil[j] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "bulk_density%02i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    soil.bulk_density[j] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "clay_in%02i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    soil.clay_in[j] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "sand_in%02i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    soil.sand_in[j] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "om%02i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    soil.om[j] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "gravel%02i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    soil.gravel[j] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "theta%02i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    soil.theta[j][1] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "theta1%02i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    wiso.dtheta[j][2] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "theta2%02i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    wiso.dtheta[j][3] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "theta3%02i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    wiso.dtheta[j][4] = (double) atof(parameter_value[i]);
	} // parameter
    } // n_soil

  // third loop through parameters, 4 times for leaf optical properties
  for (j=1; j<=4; j++)
    {
      for (i=1; i<=nparameter; i++)
	{
	  sprintf(stmp, "par_reflect_%1i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    par_reflect[j] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "par_trans_%1i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    par_trans[j] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "par_soil_refl_dry_%1i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    par_soil_refl_dry[j] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "nir_reflect_%1i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    nir_reflect[j] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "nir_trans_%1i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    nir_trans[j] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "nir_soil_refl_dry_%1i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    nir_soil_refl_dry[j] = (double) atof(parameter_value[i]);
	} // parameter
    } // leaf optical properties

  // fourth loop through parameters, betasze times for lai beta distribution
  for (j=1; j<=betasze-1; j++)
    {
      for (i=1; i<=nparameter; i++)
	{
	  sprintf(stmp, "ht_midpt%02i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    ht_midpt[j] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "lai_freq%02i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    lai_freq[j] = (double) atof(parameter_value[i]);
	  sprintf(stmp, "pai_freq%02i", j);
	  if (strcmp(parameter_name[i], stmp) == 0)
	    pai_freq[j] = (double) atof(parameter_value[i]);
	} // parameter
    } // lai beta-distribution

  // Distribute variables per layer
  if (extra_nate == 1)
    {
      for (i=1; i<=nup-1; i++)
	{
	  vc25[i] = vc25_down;
	  jm_vc[i] = jm_vc_down;
	  rd_vc[i] = rd_vc_down;
	  g0[i] = g0_down;
	  a1[i] = a1_down;
	  D0[i] = D0_down;
	  kball[i] = kball_down;
	  bprime[i] = bprime_down;
	}
      for (i=nup; i<=jtot; i++)
	{
	  vc25[i] = vc25_up;
	  jm_vc[i] = jm_vc_up;
	  rd_vc[i] = rd_vc_up;
	  g0[i] = g0_up;
	  a1[i] = a1_up;
	  D0[i] = D0_up;
	  kball[i] = kball_up;
	  bprime[i] = bprime_up;
	}
    }
  else
    {
      for (i=2; i<=jtot; i++)
	{
	  vc25[i] = vc25[1];
	  jm_vc[i] = jm_vc[1];
	  rd_vc[i] = rd_vc[1];
	  g0[i] = g0[1];
	  a1[i] = a1[1];
	  D0[i] = D0[1];
	  kball[i] = kball[1];
	  bprime[i] = bprime[1];
	}
    }

  // Calculate variables per layer
  for (i=1; i<=jtot; i++)
    {
      // calculate vcmax and jmax at optimun temperature from vcmax and jmax at 25 deg C
      jm25[i] = jm_vc[i]*vc25[i]; // calculate jmax from vcmax
      vcopt[i] = INV_TBOLTZ(vc25[i], evc, toptvc, TN0+25.);
      jmopt[i] = INV_TBOLTZ(jm25[i], ejm, toptjm, TN0+25.);
    }

}


// ----------------------------------------------------
void THROUGHFALL()
{
  /* function that calculates throughfall based on

  input.ppt : precipitation above canopy [mm]
  prof.dLAIdz: LAI per layer [m2 m-2]
  markov : clumping facotr []
  tau_water : interception efficiency
  cws_max : maximum canopy water storage per m2 LAI

  output
  prof.throughfall : throughfall [mm]  that leaves this layer
  prof.intercept : interception [mm]   that is intercepted in this layer
  prof.cws : canopy water storage [mm] that is stored in this layers leaves
  */

  int i, mc;           // counting variable for layers
  double cws_max=0.;   // maximum canopy water storage per m2 LAI [mm m-2]
  double drip=0.;      // water dripping of leaves after canopy water storage has reached max [mm]
  double intercept=0.; // interception of rain water in layer [mm]
  
  double rthrough[watersze]; // isotope ratios
  double rcws[watersze];

  cws_max = water_film_thickness*2.;   // leaves can be wet from both sides [mm m-2 LAI]

  // at canopy top: throughfall = precipitation, drip = 0
  for (mc=1; mc<=wiso.nwater+1; mc++)
    prof.throughfall[jtot1][mc] = input.ppt[mc];
  /* printf("TF01.01 %20.14f %20.14f\n", prof.throughfall[jtot+1][1], prof.throughfall[jtot+1][wiso.nwater]); */
  // calculate for each layer starting from top layer (i= jtot, typically 40)
  for (i=jtot; i>=1; i--)
    {
      // intercept per layer [mm]
      intercept = __min(tau_water*markov*prof.dLAIdz[i]*prof.throughfall[i+1][1],prof.throughfall[i+1][1]);

      for (mc=2; mc<=wiso.nwater+1; mc++)
	{
#ifdef DIAG
	  if (fabs(prof.throughfall[i+1][1]) < 1e-15 &&
	      (prof.throughfall[i+1][mc] != 0. || prof.throughfall[i+1][1] != 0.))
	    printf("\nBefore Isorat 4 tracer %i @ %ld: rare %g abundant %g\n", 
		   mc, time_var.daytime, 
		   prof.throughfall[i+1][mc], prof.throughfall[i+1][1]);
#endif
	  rthrough[mc] = ISORAT(&prof.throughfall[i+1][mc], &prof.throughfall[i+1][1],
				&wiso.lost[mc], &wiso.lost[1], mc);
	  if (wiso.lost[mc] != 0. && soil.lost0 == 0)
	    {
	      printf("\nLost 07 @ tracer %i t %ld\n",mc,time_var.daytime);
	      soil.lost0 = 1;
	    }
	}
      /* printf("TF01.02 %d %20.14f %20.14f %20.14f\n", i, intercept, rthrough[2], rthrough[wiso.nwater]); */

      // add to canopy water storage in layer (sun and shade) from previous time step
      prof.cws[i][1] += intercept;

      for (mc=2; mc<=wiso.nwater+1; mc++)
	{
	  prof.cws[i][mc] += rthrough[mc] * intercept;
#ifdef DIAG
	  if (fabs(prof.cws[i][1]) < 1e-15 &&
	      (prof.cws[i][mc] != 0. || prof.cws[i][1] != 0.))
	    printf("\nBefore Isorat 5 tracer %i @ %ld: rare %g abundant %g\n", 
		   mc, time_var.daytime, 
		   prof.cws[i][mc], prof.cws[i][1]);
#endif
	  rcws[mc] = ISORAT(&prof.cws[i][mc], &prof.cws[i][1],
			    &wiso.lost[mc], &wiso.lost[1], mc);
	  if (wiso.lost[mc] != 0. && soil.lost0 == 0)
	    {
	      printf("\nLost 08 @ tracer %i t %ld\n",mc,time_var.daytime);
	      soil.lost0 = 1;
	    }
	}
      /* printf("TF01.04 %20.14f %20.14f\n", prof.cws[i][1], prof.cws[i][wiso.nwater]); */
      /* printf("TF01.05 %20.14f %20.14f\n", rcws[2], rcws[wiso.nwater]); */

      // check if canopy water storage is at maximum
      if (prof.cws[i][1] > cws_max*prof.dLAIdz[i])
        {
          drip = prof.cws[i][1] - cws_max*prof.dLAIdz[i];
          prof.cws[i][1] = cws_max*prof.dLAIdz[i];
	  for (mc=2; mc<=wiso.nwater+1; mc++)
	    prof.cws[i][mc] = rcws[mc]*prof.cws[i][1];
        }
      else drip = 0.;
      /* printf("TF01.06 %20.14f %20.14f %20.14f\n", drip, prof.cws[i][1], prof.cws[i][wiso.nwater]); */

      // throughfall from this layer to next layer
      prof.throughfall[i][1] = prof.throughfall[i+1][1]-intercept+drip;
      for (mc=2; mc<=wiso.nwater+1; mc++)
	prof.throughfall[i][mc] = prof.throughfall[i+1][mc]-rthrough[mc]*intercept+rcws[mc]*drip;

      // Check: if this happens, we have to code some thresholds here
      for (mc=1; mc<=wiso.nwater+1; mc++)
	{
	  if (prof.throughfall[i][mc] < 0.)
	    printf("\nT<0 trac %i layer %i: %e\n", mc, i, prof.throughfall[i][mc]);
	  if (prof.cws[i][mc] < 0.)
	    printf("\nC2<0 trac %i layer %i: %e\n", mc, i, prof.cws[i][mc]);
	}
      /* printf("TF01.07 %20.14f %20.14f\n", prof.throughfall[i][1], prof.throughfall[i][wiso.nwater]); */

      // set coefficient for evaporation from wet leaf surface
      if (prof.cws[i][1] > 0.)
	prof.wet_coef_filter[i] = 1.;
      else
	prof.wet_coef_filter[i] = 0.;
    }

  return;
}


// ----------------------------------------------------
int FileCopy(const char *src, const char *dst)
{
  // Function Name: FileCopy
  // Copy the file specified by "src" to the file specified by "dst".
  // Both source and destination file names must be full path/file names,
  // such as c:\\myfile.txt a:\\myfile.txt.
  int ok=1;
  int c; /* Character read/written between files */
  FILE *ipf; /* Pointer to the I/P file. FILE is a structure  defined in <stdio.h> */
  FILE *opf;

  /* Open the file - no error checking done */
  if ((ipf = fopen(src, "rb")) == NULL) return -1;
  if ((opf = fopen(dst, "wb")) == NULL) { fclose(ipf); return -1; }

  /* Read one character at a time, checking for the End of File. EOF is defined in <stdio.h>  as -1 */
  while ((c = fgetc(ipf)) != EOF)
    fputc(c, opf); /* O/P the character */

  /* Close the files. */
  fclose(ipf);
  fclose(opf);

  return (!ok);
}


// ----------------------------------------------------
void WRITE_DUMMY(int length, ... )
{
  // writes up to 50 variables in dummy file; must be double
  // usage:
  // WRITE_DUMMY(2, vpdsoil, soil.beta);
  // WRITE_DUMMY(3, vpdsoil, soil.beta, whatever);
  double write_var[50];
  int i;
  char str[51*25];
  va_list argptr;

  if (length > 50)
    {
      puts("Too much dummy variables.");
      return;
    }

  for (i=0; i<50; i++) write_var[i] = 0.;

  va_start(argptr, length);
  for (i=0; i<length; i++) write_var[i] = va_arg(argptr, double);
  va_end(argptr);

  sprintf(str, "%ld", time_var.daytime);
  for (i=0; i<length; i++)
    {
      sprintf(str, "%s,%16.10e", str, write_var[i]);
    }
  fprintf(fptr12, "%s\n", str);

  return;
}


// ----------------------------------------------------
void SET_LEAF_PHENOLOGY()
{
  //changed: now fixed to input file, old: assume leaf out starts when soil temperature reaches 13 C at 32 cm
  time_var.leafout = leaf_out; //changed: (int)(365.*(1.5707+asin(14.-soil.Temp_ref)/soil.amplitude)/6.283);

  // assume 'leaf_full' days till full leaf
  time_var.leaffull = leaf_full; // changed: old=+30;

  // leaf shedding
  time_var.leaffall = leaf_fall;
  time_var.leaffallcomplete = leaf_fall_complete;

  return;
}


// ----------------------------------------------------
void SET_SOIL_ROOT()
{
  double r1, r2; // root fraction at top and bottom of soil layers
  double root_total; // total has to sum up to 1
  int j;

  for (j=1; j<=soil.n_soil; j++)
    {
      r1 = 1. - pow(soil.root_factor, 100.*soil.z_soil[j-1]);
      r2 = 1. - pow(soil.root_factor, 100.*soil.z_soil[j]);
      soil.root[j] = r2 - r1;
    }
  /* printf("R01 %20.14f\n", soil.root_factor); */
  /* printf("R02 %20.14f %20.14f\n", soil.z_soil[1], soil.z_soil[soil.n_soil]); */
  /* printf("R03 %20.14f %20.14f\n", soil.root[1], soil.root[soil.n_soil]); */
  
  // for Nate McDowell's juniper site give root distribution
  if (extra_nate == 1)
    {
      soil.root[1]  = 0.019927536;
      soil.root[2]  = 0.027173913;
      soil.root[3]  = 0.038043478;
      soil.root[4]  = 0.052536232;
      soil.root[5]  = 0.066847826;
      soil.root[6]  = 0.084057971;
      soil.root[7]  = 0.102717391;
      soil.root[8]  = 0.122826087;
      soil.root[9]  = 0.143115942;
      soil.root[10] = 0.342753623;
    }

  // normalise to total=1 // for comparison with v3.3.8 move above fixed root distribution
  root_total = 0.;
  for (j=1; j<=soil.n_soil; j++)
    root_total += soil.root[j];
  for (j=1; j<=soil.n_soil; j++)
    soil.root[j] /= root_total;
  /* printf("R04 %20.14f\n", root_total); */
  /* printf("R05 %20.14f %20.14f\n", soil.root[1], soil.root[soil.n_soil]); */
  
  return;
}


// ----------------------------------------------------
void SET_SOIL_TEXTURE()
{
  int j;

  for (j=1; j<=soil.n_soil; j++)
    {
      soil.clay[j] = soil.clay_in[j]*soil.clay_factor; //clay fraction
      soil.sand[j] = soil.sand_in[j]*soil.sand_factor; //sand fraction

      // watmin = 0.001;  //[m3 m-3], based on Konukcu et al. (AJSR, 2004)
      if (extra_nate == 1)
	{
	  // for Nate McDowell's juniper site
	  soil.watmin[j] = 0.1; // from measurements
	}
      else
	{
	  // estimate residual volumetric water content based on data from EPA report and selfmade regression
	  soil.watmin[j] = 0.211*soil.clay[j]/100. + 0.009*soil.sand[j]/100. + 0.028;  //[m3 m-3]
	}
    }

  // soil.qinfl_max = 5./3600.; // [mm] 5 mm h-1 is typical values from literature for loamy-clay soil
  // set maximal infiltration rate depending on soil texture at soil surface
  //   (based on FAO report, selfmade regression)
  soil.qinfl_max = __max(-26.94*soil.clay[1] + 19.27*soil.sand[1] + 13.26, 1.) / 3600.;

  return;
}


// ----------------------------------------------------
void SET_SOIL_MOISTURE()
{
  int j, mc;

  // convert delta values read from parameter file to m3 m-3
  for (j=1; j<=soil.n_soil; j++)
    {
      for (mc=2; mc<=wiso.nwater+1; mc++) // not executed if no water isotopes
	soil.theta[j][mc] = soil.theta[j][1] * INVDELTA1000_H2O(wiso.dtheta[j][mc], mc); //m3 m-3
    }

  return;
}


// ----------------------------------------------------
void SET_SOIL_TIME()
{
  soil.temperature_dt = 15.; // Time step in seconds
  soil.temperature_mtime = (long int) (time_var.time_step/soil.temperature_dt); // time steps per hour

  soil.moisture_dt = 300.; // Time step in seconds
  //soil.moisture_dt = 600.; // Time step in seconds
  soil.moisture_mtime = (long int) (time_var.time_step/soil.moisture_dt); // time steps per hour

  return;
}


// ----------------------------------------------------
void SET_SOIL_TEMP()
{
  // Assign soil temperature
  int j;

  for (j=0; j<=soil.n_soil+1; j++)
    {
      soil.T_soil[j]        = soil.T_base;
      soil.T_soil_filter[j] = soil.T_base;
    }

  return;
}


// ----------------------------------------------------
void SET_SOIL_SAXTON()
{
  // Soil characteristics after Saxton & Rawls (2006)
  int i;
  double sand[soilsze], clay[soilsze], om[soilsze];
  double theta_1500t;
  double theta_33t;
  double theta_s33t;
  double alpha, gravelw; // Matric soil density/gravel density, Weight fraction of gravel
  double psi_et;  
  double df; // used in density correction for compacted soil
  double theta_temp;  // used in density correction for compacted soil

  // Sand, Clay, OM in % in parameter file
  for (i=1; i<=soil.n_soil; i++)
    {
      sand[i] = soil.sand[i]/100.;
      clay[i] = soil.clay[i]/100.;
      om[i] = soil.om[i]/100.;
    }
  
  // Set normal density wilting point (1500 kPa), field capacity (33 kPa), excess (Sat-field)
  // air entry tension (potential), saturation moisture, density
  for (i=1; i<=soil.n_soil; i++)
    {
      theta_1500t = -0.024*sand[i] + 0.487*clay[i] + 0.006*om[i]
        + 0.005*sand[i]*om[i] - 0.013*clay[i]*om[i] + 0.068*sand[i]*clay[i]
        + 0.031;
      soil.theta_1500[i] = theta_1500t + (0.14*theta_1500t - 0.02);

      theta_33t = -0.251*sand[i] + 0.195*clay[i] + 0.011*om[i]
        + 0.006*sand[i]*om[i] - 0.027*clay[i]*om[i] + 0.452*sand[i]*clay[i]
        + 0.299;
      soil.theta_33[i] = theta_33t + (1.283*theta_33t*theta_33t - 0.374*theta_33t - 0.015);

      theta_s33t = 0.278*sand[i] + 0.034*clay[i] + 0.022*om[i]
        - 0.018*sand[i]*om[i] - 0.027*clay[i]*om[i] - 0.584*sand[i]*clay[i]
        + 0.078;
      soil.theta_s33[i] = theta_s33t + (0.636*theta_s33t - 0.107);

      psi_et = -21.674*sand[i] - 27.932*clay[i] - 81.975*soil.theta_s33[i]
        + 71.121*sand[i]*soil.theta_s33[i] + 8.294*clay[i]*soil.theta_s33[i] + 14.05*sand[i]*clay[i]
        + 27.161;
      soil.psi_e[i] = psi_et + (0.02*psi_et*psi_et - 0.113*psi_et - 0.70);

      soil.theta_s[i] = soil.theta_33[i] + soil.theta_s33[i] - 0.097*sand[i] + 0.043;

      soil.rho[i] = (1.-soil.theta_s[i]) * 2.65;
    }

  // exclude for Nate McDowell's Juniper site
  if (extra_nate == 0)
    {
      // Correct for compacted soil
      for (i=1; i<=soil.n_soil; i++)
	{
	  df = soil.bulk_density[i] / soil.rho[i];
#ifdef DIAG
	  if (df < 0.9 || df > 1.3)
	    printf("Warning: density factor not in range [0.9,1.3] of Saxton & Rawls (2006) for layer %i: %f\n"
		   "         Given density: %f. Calculated density: %f\n", i, df, soil.bulk_density[i], soil.rho[i]);
#endif
	  soil.rho[i] = soil.bulk_density[i];
	  theta_temp = 1. - soil.rho[i]/2.65;
	  soil.theta_33[i] = soil.theta_33[i] - 0.2*(soil.theta_s[i] - theta_temp);
	  soil.theta_s[i] = theta_temp;
	  soil.theta_s33[i] = soil.theta_s[i] - soil.theta_33[i];
	  // limit theta_s33 to 0.5 %v and recalc theta_33
	  soil.theta_s33[i] = __max(soil.theta_s33[i], 0.005);
	  soil.theta_33[i] = soil.theta_s[i] - soil.theta_s33[i];
	}
    }

  // Moisture parameters
  for (i=1; i<=soil.n_soil; i++)
    {
      soil.big_b[i] = (log(1500.)-log(33.)) / (log(soil.theta_33[i])-log(soil.theta_1500[i]));
      soil.big_a[i] = exp(log(33.)+soil.big_b[i]*log(soil.theta_33[i]));
      soil.lambda[i] = 1./soil.big_b[i];
      soil.k_s[i] = 1930. * pow(soil.theta_s[i]-soil.theta_33[i], 3.-soil.lambda[i]);
    }

  // Correct for gravel
  for (i=1; i<=soil.n_soil; i++)
    {
      soil.rho[i] = soil.rho[i]*(1.-soil.gravel[i]/100.) + soil.gravel[i]/100.*2.65;
      alpha = soil.rho[i]/2.65;
      gravelw = soil.gravel[i]/100. / (alpha + soil.gravel[i]/100.*(1.-alpha));
      soil.k_s[i] = soil.k_s[i] * (1.-gravelw) / (1.-gravelw*(1.-3.*alpha/2.));
    }

  // Calculate root weighted field capacity and wilting point
  soil.soil_mm_33_root = 0.; // initialize soil_mm_33_root
  soil.soil_mm_1500_root = 0.; // initialize soil_mm_1500_root
  for (i=1; i<=soil.n_soil; i++)
    {
      soil.soil_mm_33_root += soil.root[i]*soil.theta_33[i]*1000.*(soil.z_soil[i]-soil.z_soil[i-1])*(1.-soil.gravel[i]/100.); // soil.soil_mm_33_root is used for soil water stress
      soil.soil_mm_1500_root += soil.root[i]*soil.theta_1500[i]*1000.*(soil.z_soil[i]-soil.z_soil[i-1])*(1.-soil.gravel[i]/100.); // soil.soil_mm_1500_root is used for soil water stress
    }


  return;
}


// ----------------------------------------------------
void SET_SOIL_CLAPP()
{
  // Soil characteristics after Clapp & Hornberger (1978)
  int i;
  double sand[soilsze], clay[soilsze];
  double alpha, gravelw; // Matric soil density/gravel density, Weight fraction of gravel

  // Sand, Clay, OM in % in parameter file
  for (i=1; i<=soil.n_soil; i++)
    {
      sand[i] = soil.sand[i];
      clay[i] = soil.clay[i];
    }

  for (i=1; i<=soil.n_soil; i++)
    {
      soil.theta_s[i] = 0.489 - 0.00126*soil.sand[i];
      soil.k_s[i] = 0.0070556 * ( pow( 10.,(-0.884+0.0153*soil.sand[i]) ) );
      soil.psi_e[i] = -10. * ( pow(10.,(1.88-0.0131*soil.sand[i])) );
      soil.big_b[i] = 2.91 + 0.159*soil.clay[i];
      soil.theta_1500[i] = soil.theta_s[i] * pow(-316230./soil.psi_e[i],-1./soil.big_b[i]);
      soil.theta_33[i] = soil.theta_s[i] * pow(-158490./soil.psi_e[i],-1./soil.big_b[i]);
      soil.theta_s33[i] = soil.theta_s[i] - soil.theta_33[i];
      // limit theta_s33 to 0.5 %v and recalc theta_33
      soil.theta_s33[i] = __max(soil.theta_s33[i], 0.005);
      soil.theta_33[i] = soil.theta_s[i] - soil.theta_s33[i];
      if (extra_nate == 1)
	{ // Calc density for Nate McDowell's Juniper site
	  soil.rho[i] = (1.-soil.theta_s[i]) * 2.65;
	}
      else
	soil.rho[i] = soil.bulk_density[i];
    }

  // Correct for gravel
  for (i=1; i<=soil.n_soil; i++)
    {
      soil.rho[i] = soil.rho[i]*(1.-soil.gravel[i]/100.) + soil.gravel[i]/100.*2.65;
      alpha = soil.rho[i]/2.65;
      gravelw = soil.gravel[i]/100. / (alpha + soil.gravel[i]/100.*(1.-alpha));
      soil.k_s[i] = soil.k_s[i] * (1.-gravelw) / (1.-gravelw*(1.-3.*alpha/2.));
    }

  // Calculate root weighted field capacity and wilting point
  soil.soil_mm_33_root = 0.; // initialize soil_mm_33_root
  soil.soil_mm_1500_root = 0.; // initialize soil_mm_1500_root
  for (i=1; i<=soil.n_soil; i++)
    {
      soil.soil_mm_33_root += soil.root[i]*soil.theta_33[i]*1000.*(soil.z_soil[i]-soil.z_soil[i-1])*(1.-soil.gravel[i]/100.); // soil.soil_mm_33_root is used for soil water stress
      soil.soil_mm_1500_root += soil.root[i]*soil.theta_1500[i]*1000.*(soil.z_soil[i]-soil.z_soil[i-1])*(1.-soil.gravel[i]/100.); // soil.soil_mm_1500_root is used for soil water stress
    }


  return;
}


// ----------------------------------------------------
double SOIL_RESISTANCE()
{
  double c = 13.515; // soil resistance of Passerat (1986), from Mahfouf & Noilhan (1991)
  double rsmax = 3.8113e4;
  double dw, dw_0 = 0.229e-4; // vapour diffusion in air @ TN0;
  double rs, rl;

  rs = rsmax * exp(-c * (soil.theta[1][1]-soil.watmin[1]) / (soil.theta_33[1]-soil.watmin[1]));

  // Camillo and Gurney model for soil resistance
  // Rsoil= 4104 (ws-wg)-805, ws=.395, wg=0
  // ws = max soil moisture = 0.55 at the Hainich site, wg=0
  if (soil.camillo == 1) rs = 1452.2;

  if (soil.z_litter < 1e-6)
    rl = 0.;
  else
    {
      dw = dw_0 * pow(soil.T_l_filter/TN0, 1.75); // from Kondo et al. (1990)
      rl = soil.z_litter/(dw*pow(__max(soil.theta_ls-soil.theta_l[1], 1e-12), 5./3.)); // with Milly model (5/3) after Saravanapavan et al. (2000).
    }

  return (rs+rl);
}


// ----------------------------------------------------
double SOIL_RH()
{
  double s, smp, rh;

  if (soil.camillo == 1)
    rh = 1.;
  else
    {
      if (soil.saxton == 1)
        {
          if (soil.theta[1][1] <= soil.theta_33[1]) // swp in kPa
            smp = soil.big_a[1] * pow(soil.theta[1][1], -soil.big_b[1]);
          else
            smp = 33. - (33.-soil.psi_e[1])*(__min(soil.theta[1][1], soil.theta_s[1])-soil.theta_33[1])
              /(soil.theta_s[1]-soil.theta_33[1]);
          smp = -smp*1000./Gravity; // kPa -> mm
        }
      else
        {
          s = __max(soil.theta[1][1]/soil.theta_s[1], 0.05); // % saturation
          smp = soil.psi_e[1] * pow(s, -soil.big_b[1]); // soil water potential swp (mm water)
        }
      rh = exp(-Gravity*smp/1000./(Rw*(soil.T_soil_filter[1]+TN0)));
    }

  return __max(__min(rh, 1.), 0.);
}


// ----------------------------------------------------
void SET_LITTER_TEXTURE()
{
  double rho_0 = 80.; // dry density of std organic matter [kg m-3]
                      // average of table 2 of Redding et al. (Soil Science Society of America Journal, 2005)
  soil.theta_ls  = 1. - soil.rho_l/rho_0; // same as Saxton & Rawls (2006), i.e. theta_s = porosity
  soil.theta_l33 = 0.25; // fit of rmax*(1-theta/theta_field_capacity)^n to data of Schaap & Bouten (1997)
  soil.n_l       = 2.5;        // dito, see fit_litter_resistance.pro

  return;
}


// ----------------------------------------------------
void SET_LITTER_TEMP()
{
  // Assign litter temperature
  soil.T_l        = soil.T_base+TN0;
  soil.T_l_filter = soil.T_base+TN0;

  return;
}


// ----------------------------------------------------
void SET_LITTER_MOISTURE()
{
  int mc;

  if (soil.z_litter < 1e-6)
    for (mc=1; mc<=wiso.nwater+1; mc++)
      soil.theta_l[mc] = 0.; //m3 m-3
  else
    for (mc=1; mc<=wiso.nwater+1; mc++)
      soil.theta_l[mc] = soil.theta[1][mc]; //m3 m-3

  return;
}


// ----------------------------------------------------
double LITTER_RESISTANCE()
{
  double dw, dw_0 = 0.229e-4; // vapour diffusion in air @ TN0;
  double rlmax, rl;

  dw = dw_0 * pow(soil.T_l_filter/TN0, 1.75); // from Kondo et al. (1990)  
  if (extra_nate == 1)
    {
      // for Nate McDowell's juniper site, fix maximum litter resistance
      rlmax = 5050;    //soil.z_litter/(dw* pow(soil.theta_ls, 2./3.)); // after Schaap & Bouten (1997)
      rlmax = 4200;
    }
  else
    {
      rlmax = soil.z_litter/(dw* pow(soil.theta_ls, 2./3.)); // after Schaap & Bouten (1997)
      // with Milly model (2/3) after Saravanapavan et al. (2000).
    }

  rl = 0.;
  if (soil.theta_l[1] < soil.theta_l33)
    rl = rlmax*pow(1.-soil.theta_l[1]/soil.theta_l33, soil.n_l);

  return rl;
}


// ----------------------------------------------------
void LITTER_RAIN()
{
  double dw_l[watersze];
  int mc;

  for (mc=1; mc<=wiso.nwater+1; mc++)
    flux.soilinfl[mc] = prof.throughfall[1][mc]/time_var.time_step;  // convert from mm to mm/s
  
  // new litter water status, mass balance
  dw_l[1] = flux.soilinfl[1] / 1000. / soil.z_litter; // mm/s
  soil.theta_l[1] = soil.theta_l[1] + dw_l[1]*time_var.time_step; // mm
  for (mc=2; mc<=wiso.nwater+1; mc++)
    {
      dw_l[mc] = flux.soilinfl[mc] / 1000. / soil.z_litter;
      soil.theta_l[mc] = soil.theta_l[mc] + dw_l[mc]*time_var.time_step;
    }

  return;
}


// ----------------------------------------------------
void LITTER_H2O()
{
  double rsoil, dw_l[watersze];
  int mc;
  double tmp1, tmp2, rqinfl;

  if (soil.z_litter < 1e-6)
    {
      for (mc=1; mc<=wiso.nwater+1; mc++)
	flux.soilinfl[mc] = prof.throughfall[1][mc]/time_var.time_step;  // convert from mm to mm/s

      // if z_litter = 0 then theta_l is something like a rain buffer or a puddle or similar [mm]
      tmp1 = flux.soilinfl[1]+soil.theta_l[1]/time_var.time_step;
      soil.qinfl[1] = __min(tmp1, soil.qinfl_max); // [kg/m2/s = mm/s]

      for (mc=2; mc<=wiso.nwater+1; mc++)
	{
	  tmp2 = flux.soilinfl[mc]+soil.theta_l[mc]/time_var.time_step;
#ifdef DIAG
	  if (fabs(tmp1) < 1e-15 &&
	      (tmp1 != 0. || tmp2 != 0.))
	    printf("\nBefore Isorat 6 tracer %i @ %ld: rare %g abundant %g\n", 
		   mc, time_var.daytime, 
		   tmp2, tmp1);
#endif
	  rqinfl = ISORAT(&tmp2, &tmp1, &wiso.lost[mc], &wiso.lost[1], mc);
	  if (wiso.lost[mc] != 0. && soil.lost0 == 0)
	    {
	      printf("\nLost 09 @ tracer %i t %ld\n",mc,time_var.daytime);
	      soil.lost0 = 1;
	    }
	  soil.qinfl[mc] = soil.qinfl[1] * rqinfl;
	}

      soil.theta_l[1] = soil.theta_l[1]
	+ flux.soilinfl[1]*time_var.time_step
	- soil.qinfl[1]*time_var.time_step;

      for (mc=2; mc<=wiso.nwater+1; mc++)
	{
	  soil.theta_l[mc] = soil.theta_l[mc]
	    + flux.soilinfl[mc]*time_var.time_step
	    - soil.qinfl[mc]*time_var.time_step;
	}

      // Check: if this happens, we have to code some thresholds here
      for (mc=1; mc<=wiso.nwater+1; mc++)
	{
	  if (soil.theta_l[mc] < (-isotolerance))
	    printf("\nL1<0 trac %i: %e\n", mc, soil.theta_l[mc]);
	  if (soil.theta_l[mc] < isotolerance)
	    {
	      soil.qinfl[mc] += soil.theta_l[mc];
	      soil.theta_l[mc] = 0.;
	    }
	}
    }
  else
    {
      //  Litter water drainage into soil, from Ogee & Brunet (2002)
      soil.qinfl[1] = __min(__max(0., 
					 1000.*soil.z_litter
					 *(soil.theta_l[1]-soil.theta_l33)
					 /time_var.time_step),
				   soil.qinfl_max); //  kg/m2/s = mm/s

      for (mc=2; mc<=wiso.nwater+1; mc++)
        {
#ifdef DIAG
	  if (fabs(soil.theta_l[1]) < 1e-15 &&
	      (soil.theta_l[mc] != 0. || soil.theta_l[1] != 0.))
	    printf("\nBefore Isorat 7 tracer %i @ %ld: rare %g abundant %g\n", 
		   mc, time_var.daytime, 
		   soil.theta_l[mc], soil.theta_l[1]);
#endif
          rsoil = ISORAT(&soil.theta_l[mc], &soil.theta_l[1],
			 &wiso.lost[mc], &wiso.lost[1], mc);
	  if (wiso.lost[mc] != 0. && soil.lost0 == 0)
	    {
	      printf("\nLost 10 @ tracer %i t %ld\n",mc,time_var.daytime);
	      soil.lost0 = 1;
	    }
          soil.qinfl[mc] = rsoil * soil.qinfl[1];
        }

      // new litter water status, mass balance
      dw_l[1] = (0. - flux.litterevap[1] - soil.qinfl[1]) / 1000. / soil.z_litter;
      if ((soil.theta_l[1]+dw_l[1]*time_var.time_step) < 0.)
	{
	  printf("%ld L2<<0 %i %e %e %e %e\n", time_var.daytime, 1
		 , flux.litterevap[1]/LAMBDA(soil.T_l)/1000./soil.z_litter*time_var.time_step
		 , soil.qinfl[1]/1000./soil.z_litter*time_var.time_step
		 , soil.theta_l[1]
		 , dw_l[1]*time_var.time_step);
	}
      soil.theta_l[1] = soil.theta_l[1] + dw_l[1]*time_var.time_step;
      for (mc=2; mc<=wiso.nwater+1; mc++)
        {
          dw_l[mc] = (0. - flux.litterevap[mc] - soil.qinfl[mc])
	    / 1000. / soil.z_litter;
	  if ((soil.theta_l[mc]+dw_l[mc]*time_var.time_step) < 0.)
	    {
	      printf("%ld L2<<0 %i l %e q %e t %e d %e l1 %e q1 %e t1 %e\n", time_var.daytime, mc
		     , flux.litterevap[mc]/LAMBDA(soil.T_l)/1000./soil.z_litter*time_var.time_step
		     , soil.qinfl[mc]/1000./soil.z_litter*time_var.time_step
		     , soil.theta_l[mc]
		     , dw_l[mc]*time_var.time_step
		     , flux.litterevap[1]/LAMBDA(soil.T_l)/1000./soil.z_litter*time_var.time_step
		     , soil.qinfl[1]/1000./soil.z_litter*time_var.time_step
		     , soil.theta_l[1]
		     );
	    }
          soil.theta_l[mc] = soil.theta_l[mc] + dw_l[mc]*time_var.time_step;
        }

      // Check: if this happens, we have to code some thresholds here
      for (mc=1; mc<=wiso.nwater+1; mc++)
	if (soil.theta_l[mc] < 0.)
	  printf("\nL2<0 trac %i: %e\n", mc, soil.theta_l[mc]);
    }

  return;
}


// ----------------------------------------------------
double LITTER_RH()
{
  double psi_l1 = 35.3; // Ogee & Brunet (2002)
  double bl = 2.42;
  double rho_w = 1000.; // density of water
  double psi_l, rh;

  if (soil.theta_l[1] != 0.)
    psi_l = psi_l1 * pow(rho_w*soil.theta_l[1]/soil.rho_l, -bl);
  else
    psi_l = 0.;
  rh = exp(-Gravity*psi_l/(Rw*soil.T_l_filter));
  
  return __max(__min(rh, 1.), 0.);
}


// ----------------------------------------------------
double ISORAT(double *qi, double *q, double *lostqi, double *lostq, int species)
{
  // returns isotope ratio qi/q
  // checks for underflow. If abs(q)<1e-15 -> lostqi=lostqi+qi and qi=q=0
  // double precision allows for about 15 digits of precision -> 1e-15
  //MC20190627 double lower_bound = 1e-15;
  double lower_bound = 2.2204460492503131E-016; // epsilon(1.) of Fortran
  double tiny = 2.2250738585072014E-308;        // tiny(1.) of Fortran
  double ratio;

  if (fabs(*q) >= lower_bound)
    ratio = *qi / *q;
  else
    {
      //MC20190627 if (*q != 0. || *qi != 0.)
      if (fabs(*q) > tiny || fabs(*qi) > tiny)
	{
	  printf("\nIsorat tracer %i @ %ld: rare %g abundant %g\n", species, time_var.daytime, *qi, *q);
	  //soil.lost0 = 1;
	}
      *lostqi += *qi;
      *qi = 0.;
      *lostq += *q;
      *q = 0.;
      ratio = 0.;
    }

  return ratio;
}


// ----------------------------------------------------
int TRIDIA ()
{

  /* ------------------------ code history ---------------------------
   * source file: tridia.F
   * purpose: solve tridiagonal system of equations
   * date last revised: March 1996 - lsm version 1
   * author: Gordon Bonan
   * standardized: J. Truesdale, Feb. 1996
   * reviewed: G. Bonan, Feb. 1996
   * -----------------------------------------------------------------

   * ------------------------ notes ----------------------------------
   * solve for u given the set of equations f * u = r, where u is a
   * vector of length n, r is a vector of length n, and f is an n x n
   * tridiagonal matrix defined by the vectors a, b, c [each of length n].
   * a(1) and c(n) are undefined and are not referenced by the subroutine.

   * |b(1) c(1) 0 ... | |u(1 )| |r(1 )|
   * |a(2) b(2) c(2) ... | |u(2 )| |r(2 )|
   * | ... | * | ... | = | ... |
   * | ... a(n-1) b(n-1) c(n-1)| |u(n-1)| |r(n-1)|
   * | ... 0 a(n ) b(n )| |u(n )| |r(n )|

   * -----------------------------------------------------------------*/

  int ok=1;

  //* ------------------------ input/output variables -----------------
  double a[soilsze]; // input vector
  double b[soilsze]; // input vector
  double c[soilsze]; // input vector
  double r[soilsze][watersze]; // input vector
  double u[soilsze]; // solution vector
  //* -----------------------------------------------------------------

  //* ------------------------ local variables ------------------------
  double gam[soilsze];
  double bet;
  int j, mc; // loop index
  //* -----------------------------------------------------------------

  for (mc=1; mc<=wiso.nwater+1; mc++)
    {
      for (j=1; j<=soil.n_soil; j++) //old just <
        {
          r[j][mc] = soil.r[j][mc];
          a[j] = soil.a[j];
          b[j] = soil.b[j];
          c[j] = soil.c[j];
        }

      bet = b[1];
      u[1] = r[1][mc] / bet;

      for (j=2; j<=soil.n_soil; j++) //old just <
        {
          gam[j] = c[j-1] / bet;
          bet = b[j] - a[j] * gam[j];
          u[j] = (r[j][mc] - a[j]*u[j-1]) / bet;
        }

      for (j=soil.n_soil-1; j>=1; j--) //old just <
        {
          u[j] = u[j] - gam[j+1] * u[j+1];
        }

      for (j=1; j<=soil.n_soil; j++)
        {
          soil.dwat[j][mc] = u[j];
        }
    }
  return (!ok);
}


// ----------------------------------------------------
int IUP(int ihere, double flx)
{
  if (flx <= 0.)
    return ihere;
  else
    return ihere+1;
}


// ----------------------------------------------------
int SOIL_H2O()
{

  int ok=1;
  /* ------------------------ code history ------------------------------
   * source file: soih2o.F
   * purpose: soil hydrology and sub-surface drainage
   * date last revised: March 1996 - lsm version 1
   * author: Gordon Bonan
   * standardized: J. Truesdale, Feb. 1996
   * reviewed: G. Bonan, Feb. 1996
   * --------------------------------------------------------------------

   * ------------------------ notes -------------------------------------
   * use tridiagonal system of equations to solve one-dimensional water balance:

   * d wat
   * dz ----- = -qi + qo - s
   * dt

   * with q = -k d(psi+z)/dz = -k (d psi/dz + 1) and with s=0
   * this is the Richards equation for vertical water flow

   * d wat d d wat d psi
   * ----- = -- [ k(----- ----- + 1) ]
   * dt dz dz d wat

   * where: wat = volume of water per volume of soil (mm**3/mm**3)
   * psi = soil matrix potential (mm)
   * dt = time step (s)
   * z = depth (mm)
   * dz = thickness (mm)
   * qi = inflow at top (mm h2o /s) (+ = up, - = down)
   * qo = outflow at bottom (mm h2o /s) (+ = up, - = down)
   * s = source/sink flux (mm h2o /s) (+ = loss, - = gain)
   * k = hydraulic conductivity (mm h2o /s)

   * solution: linearize k and psi about d wat and use tridiagonal system
   * of equations to solve for d wat, where for layer i
   * r_i = a_i [d wat_i-1] + b_i [d wat_i] + c_i [d wat_i+1]

   * the solution conserves water as:
   * [h2osoi(1)*dzsoi(1)*1000 + ... + h2osoi(nsl)*dzsoi(nsl)*1000] n+1 =
   * [h2osoi(1)*dzsoi(1)*1000 + ... + h2osoi(nsl)*dzsoi(nsl)*1000] n +
   * (qinfl - qseva - qtran - qdrai)*dtlsm

   * --------------------------------------------------------------------
   */
  //* ------------------------ input/output variables --------------------
  //* input
  double qtran[watersze]; //transpiration water flux (mm h2o/s)
  double qinfl[watersze]; //infiltration rate (mm h2o/s)
  double qseva[watersze]; //ground surface evaporation rate (mm h2o/s)
  double dzsoi[soilsze]; //soil layer thickness (m)

  double root[soilsze]; //relative root abundance (0 to 1)

  //* input/output
  double h2osoi[soilsze][watersze]; //volumetric soil water content
  double lost_h2osoi[watersze]; // water conservation check variable

  //* output
  double qdrai[watersze]; //sub-surface runoff (mm h2o/s)

  //* --------------------------------------------------------------------

  //* ------------------------ local variables ---------------------------
  int j, mc; //do loop/array indices; mc for water isotopes
  int iter; //loop

  double r[soilsze][watersze]; //solution matrix
  double a[soilsze]; //"a" vector for tridiagonal matrix
  double b[soilsze]; //"b" vector for tridiagonal matrix
  double c[soilsze]; //"c" vector for tridiagonal matrix
  double dwat[soilsze][watersze]; //change in soil water
  double smp[soilsze]; //soil matrix potential (mm)
  double hk[soilsze]; //hydraulic conductivity (mm h2o/s)
  double hk2[soilsze]; //hk**2
  double dsmpdw[soilsze]; //d(smp)/d(h2osoi)
  double dhkdw[soilsze]; //d(hk)/d(h2osoi)
  double s; // relative soil water content in Clapp & Hornberger
  double hydcon; // temp hydraulic conductivity

  double qin; //flux of water into soil layer (mm h2o/s)
  double qout; //flux of water out of soil layer (mm h2o/s)
  double num; //used in calculating qi, qout
  double den; //used in calculating qi, qout
  double den2; //den**2 used in calculating qi, qout
  double dqidw1; //d(qin)/d(h2osoi(i-1))
  double dqidw2; //d(qin)/d(h2osoi(i))
  double dqodw1; //d(qout)/d(h2osoi(i))
  double dqodw2; //d(qout)/d(h2osoi(i+1))
  double xs[watersze]; //soil water > sat or < some minimum (mm h2o)

  double thisxs[watersze];
  double rxs, rdrai;
  double hkdrai;
  double rsoil[soilsze][watersze]; // ratio of levels

  double sm, sm50;
  double xstmp1[soilsze][watersze];
  double xstmp2[soilsze][watersze];

  //* --------------------------------------------------------------------

  // Assign local variables

  for (mc=1; mc<=wiso.nwater+1; mc++) // not executed if no water isotopes
    {
      soil.qtran[mc] = qtran[mc] = flux.c_transpiration[mc]; //transpiration water flux (kg H2O m-2 s-1) = mm h2o/s
      qinfl[mc] = soil.qinfl[mc];
      soil.qseva[mc] = qseva[mc] = flux.soilevap[mc]; //ground surface evaporation rate (mm h2o/s) = (kg H2O m-2 s-1)
      soil.qdrai[mc] = qdrai[mc] = 0.;
    }

  for (j=1; j<=soil.n_soil; j++)
    {
      dzsoi[j] = (soil.z_soil[j]-soil.z_soil[j-1])*(1.-soil.gravel[j]/100.); // set effective soil layer thickness
      root[j] = soil.root[j]; //relative root abundance (0 to 1)
      for (mc=1; mc<=wiso.nwater+1; mc++)
        h2osoi[j][mc] = soil.theta[j][mc]; //volumetric soil water content, get soil water content from last time step
    }
  for (mc=1; mc<=wiso.nwater+1; mc++)
    lost_h2osoi[mc] = wiso.lost[mc];

  for (iter=1; iter<=soil.moisture_mtime; iter++)
    {

      for (j=1; j<=soil.n_soil; j++)
        for (mc=2; mc<=wiso.nwater+1; mc++)
	  {
#ifdef DIAG
	    if (fabs(h2osoi[j][1]) < 1e-15 &&
		(h2osoi[j][mc] != 0. || h2osoi[j][1] != 0.))
	      printf("\nBefore Isorat 8 tracer %i @ %ld: rare %g abundant %g\n", 
		     mc, time_var.daytime, 
		     h2osoi[j][mc], h2osoi[j][1]);
#endif
	    rsoil[j][mc] = ISORAT(&h2osoi[j][mc], &h2osoi[j][1], &lost_h2osoi[mc], &lost_h2osoi[1], mc);
	    if (lost_h2osoi[mc] != 0. && soil.lost0 == 0)
	      {
		printf("\nLost 11.1 @ tracer %i t %ld\n", mc, time_var.daytime);
		soil.lost0 = 1;
	      }
	  }

      /* evaluate hydraulic conductivity, soil matrix potential,
       * d(smp)/d(h2osoi), and d(hk)/d(h2osoi). limit s >= 0.05
       * when evaluating these terms. this helps prevent numerical
       * problems for very small h2osoi. also limit hk >= some very
       * small number for same reason.*/

      if (soil.saxton == 1)
        {
          for (j=1; j<=soil.n_soil; j++)
            {
              if (h2osoi[j][1] <= soil.theta_33[j]) // swp in kPa
                smp[j] = soil.big_a[j] * pow(__max(h2osoi[j][1],soil.watmin[j]), -soil.big_b[j]);
              else
                smp[j] = 33. - (33.-soil.psi_e[j])
		  * (__min(__max(h2osoi[j][1],soil.watmin[j]), soil.theta_s[j])
		     -soil.theta_33[j])
		  / (soil.theta_s[j]-soil.theta_33[j]);
              soil.swp[j] = smp[j]; //kPa

              if (h2osoi[j][1] <= soil.theta_33[j]) // derivative of swp to sm
		// same as Clapp & Hornberger
                dsmpdw[j] = -soil.big_b[j]*smp[j]/__max(h2osoi[j][1],soil.watmin[j]);
              else
                dsmpdw[j] = - (33.-soil.psi_e[j])/(soil.theta_s[j]-soil.theta_33[j]);

              soil.swp_mm[j] = smp[j] = -smp[j]*1000./Gravity; // kPa -> mm
              dsmpdw[j] = -dsmpdw[j]*1000./Gravity; // kPa -> mm

              soil.k_theta[j] = hk[j] 
		= soil.k_s[j] 
		* pow(__max(h2osoi[j][1],soil.watmin[j])/soil.theta_s[j],
		      3.+2./soil.lambda[j]); // mm h-1
	      // derivative of hk on sm
              dhkdw[j] = (3.+2./soil.lambda[j])*hk[j]/__max(h2osoi[j][1],soil.watmin[j]);
              soil.k_theta[j] = hk[j] = hk[j] / 3600.; // mm s-1
              dhkdw[j] = dhkdw[j] / 3600.;
              hk2[j] = hk[j]*hk[j]; // hk^2
            }
        }
      else
        {
          for (j=1; j<=soil.n_soil; j++) // Clapp & Hornberger
            {
              s = __max(h2osoi[j][1]/soil.theta_s[j],0.05); // % de saturation
              smp[j] = soil.psi_e[j] * pow(s, -soil.big_b[j]); // soil water potential swp (mm water)
              soil.swp_mm[j] = smp[j]; // mm
              soil.swp[j] = -smp[j]*Gravity/1000.; // kPa
	      // derivative of swp on sm
              dsmpdw[j] = -soil.big_b[j]*smp[j]/__max(h2osoi[j][1],soil.watmin[j]);
              hydcon = soil.k_s[j] * pow(s, 2.*soil.big_b[j]+3.); // hydraulic conductivity (mm h2o/s)
              hk[j] = __max(hydcon, 1.e-10); // hydraulic conductivity hk is not < 0.0
	      // derivative of hk on sm
              dhkdw[j] = hydcon*(2.*soil.big_b[j]+3.)/__max(h2osoi[j][1],soil.watmin[j]);
              hk2[j] = hk[j]*hk[j]; // hk^2
            }
        }

      //* set up r, a, b, and c vectors for tridiagonal solution

      //* node j=1

      j = 1;
      num = -2.*(smp[j]-smp[j+1]) - 1000.*(dzsoi[j]+dzsoi[j+1]);
      den = 1000.*(dzsoi[j]/hk[j] + dzsoi[j+1]/hk[j+1]);
      den2 = den*den;
      qout = num/den;

      dqodw1 = (-2.*den*dsmpdw[j ] + num*1000.*dzsoi[j ]/hk2[j ]*dhkdw[j ])/ den2;
      dqodw2 = ( 2.*den*dsmpdw[j+1] + num*1000.*dzsoi[j+1]/hk2[j+1]*dhkdw[j+1])/ den2;

      mc = 1;
      r[j][mc] = (qseva[mc]+qtran[mc]*root[j]) - qinfl[mc] - qout;
      for (mc=2; mc<=wiso.nwater+1; mc++)
        {
#ifdef DIAG
	  if (fabs(h2osoi[j][1]) < 1e-15 &&
	      (h2osoi[j][mc] != 0. || h2osoi[j][1] != 0.))
	    printf("\nBefore Isorat 9 tracer %i @ %ld: rare %g abundant %g\n", 
		   mc, time_var.daytime, 
		   h2osoi[j][mc], h2osoi[j][1]);
#endif
          r[j][mc] = qseva[mc] +
            ISORAT(&h2osoi[j][mc], &h2osoi[j][1], &lost_h2osoi[mc], &lost_h2osoi[1], mc) *
            qtran[1]*root[j] -
            qinfl[mc] - qout*rsoil[IUP(j,qout)][mc];
	  if (lost_h2osoi[mc] != 0. && soil.lost0 == 0)
	    {
	      printf("\nLost 11.2 @ tracer %i t %ld\n", mc, time_var.daytime);
	      soil.lost0 = 1;
	    }
        }
      a[j] = 0.;
      b[j] = dqodw1 - 1000.*dzsoi[j]/soil.moisture_dt;
      c[j] = dqodw2;

      //* node j=n_soil

      j = soil.n_soil;
      num = -2.*(smp[j-1]-smp[j]) -1000.*(dzsoi[j-1]+dzsoi[j]);
      den = 1000.*(dzsoi[j-1]/hk[j-1] + dzsoi[j]/hk[j]);
      den2 = den*den;
      qin = num/den;

      dqidw1 = (-2.*den*dsmpdw[j-1] +num*1000.*dzsoi[j-1]/hk2[j-1]*dhkdw[j-1])/ den2;
      dqidw2 = ( 2.*den*dsmpdw[j ] +num*1000.*dzsoi[j ]/hk2[j ]*dhkdw[j ])/ den2;

      qout = -hk[j];
      dqodw1 = -dhkdw[j];

      mc = 1;
      r[j][mc] = qtran[1]*root[j] + qin - qout;
      for (mc=2; mc<=wiso.nwater+1; mc++)
        {
#ifdef DIAG
	  if (fabs(h2osoi[j][1]) < 1e-15 &&
	      (h2osoi[j][mc] != 0. || h2osoi[j][1] != 0.))
	    printf("\nBefore Isorat 10 tracer %i @ %ld: rare %g abundant %g\n", 
		   mc, time_var.daytime, 
		   h2osoi[j][mc], h2osoi[j][1]);
#endif
          r[j][mc] = ISORAT(&h2osoi[j][mc], &h2osoi[j][1], &lost_h2osoi[mc], &lost_h2osoi[1], mc) *
            qtran[1]*root[j] +
            qin*rsoil[IUP(j-1,qin)][mc] -
            qout*rsoil[soil.n_soil][mc];
	  if (lost_h2osoi[mc] != 0. && soil.lost0 == 0)
	    {
	      printf("\nLost 11.3 @ tracer %i t %ld\n", mc, time_var.daytime);
	      soil.lost0 = 1;
	    }
        }
      a[j] = -dqidw1;
      b[j] = dqodw1 - dqidw2 - 1000.*dzsoi[j]/soil.moisture_dt;
      c[j] = 0.;

      //* nodes j=2 to j=n_soil-1

      for (j=2; j<=soil.n_soil-1; j++)
        {
          num = -2.*(smp[j-1]-smp[j]) -1000.*(dzsoi[j-1]+dzsoi[j]);
          den = 1000.*(dzsoi[j-1]/hk[j-1] + dzsoi[j]/hk[j]);
          den2 = den*den;
          qin = num/den;
          dqidw1 = (-2.*den*dsmpdw[j-1] +num*1000.*dzsoi[j-1]/hk2[j-1]*dhkdw[j-1])/ den2;
          dqidw2 = ( 2.*den*dsmpdw[j ] +num*1000.*dzsoi[j ]/hk2[j ]*dhkdw[j ])/ den2;

          num = -2.*(smp[j]-smp[j+1]) -1000.*(dzsoi[j]+dzsoi[j+1]);
          den = 1000.*(dzsoi[j]/hk[j] + dzsoi[j+1]/hk[j+1]);
          den2 = den*den;
          qout = num/den;
          dqodw1 = (-2.*den*dsmpdw[j ] +num*1000.*dzsoi[j ]/hk2[j ]*dhkdw[j ])/ den2;
          dqodw2 = ( 2.*den*dsmpdw[j+1] +num*1000.*dzsoi[j+1]/hk2[j+1]*dhkdw[j+1])/ den2;

          mc = 1;
          r[j][mc] = qtran[mc] * root[j] + qin - qout;
          for (mc=2; mc<=wiso.nwater+1; mc++)
            {
#ifdef DIAG
	      if (fabs(h2osoi[j][1]) < 1e-15 &&
		  (h2osoi[j][mc] != 0. || h2osoi[j][1] != 0.))
		printf("\nBefore Isorat 11 tracer %i @ %ld: rare %g abundant %g\n", 
		       mc, time_var.daytime, 
		       h2osoi[j][mc], h2osoi[j][1]);
#endif
              r[j][mc] = ISORAT(&h2osoi[j][mc], &h2osoi[j][1], &lost_h2osoi[mc], &lost_h2osoi[1], mc) *
                qtran[1]*root[j] +
                qin*rsoil[IUP(j-1,qin)][mc] -
                qout*rsoil[IUP(j,qout)][mc];
	      if (lost_h2osoi[mc] != 0. && soil.lost0 == 0)
		{
		  printf("\nLost 11.4 @ tracer %i t %ld\n", mc, time_var.daytime);
		  soil.lost0 = 1;
		}
            }
          a[j] = -dqidw1;
          b[j] = dqodw1 - dqidw2 - 1000.*dzsoi[j]/soil.moisture_dt;
          c[j] = dqodw2;
        }

      //* solve for dwat: a, b, c, r, dwat

      for (j=1; j<=soil.n_soil; j++)
        {
          soil.a[j] = a[j];
          soil.b[j] = b[j];
          soil.c[j] = c[j];
          for (mc=1; mc<=wiso.nwater+1; mc++)
            {
              soil.r[j][mc] = r[j][mc];
            }
        }

      TRIDIA();

      /* could now update h2osoi = h2osoi + dwat except for one problem:
       * need to make sure h2osoi <= watsat. if not, what to do with
       * excess water? 
       * We percolate any excess water down the soil. This is different from ISOLSM.
       * Any excess water in the last soil layer is sub-surface runoff
       */
      for (j=1; j<=soil.n_soil; j++)
        for (mc=1; mc<=wiso.nwater+1; mc++)
	  {
	    dwat[j][mc] = soil.dwat[j][mc];
            h2osoi[j][mc] += dwat[j][mc];
	  }

      // from 1 to n_soil-1
      for (j=1; j<=soil.n_soil-1; j++)
        {
          thisxs[1] = __max(h2osoi[j][1]-soil.theta_s[j], 0.) * 1000. * dzsoi[j];
	  xstmp1[j][1] = thisxs[1];
	  if (thisxs[1] != 0)
	    {
	      for (mc=2; mc<=wiso.nwater+1; mc++)
		{
#ifdef DIAG
		  if (fabs(h2osoi[j][1]) < 1e-15 &&
		      (h2osoi[j][mc] != 0. || h2osoi[j][1] != 0.))
		    printf("\nBefore Isorat 12 tracer %i @ %ld: rare %g abundant %g\n", 
			   mc, time_var.daytime, 
			   h2osoi[j][mc], h2osoi[j][1]);
#endif
		  thisxs[mc] = thisxs[1] *
		    ISORAT(&h2osoi[j][mc], &h2osoi[j][1], &lost_h2osoi[mc], &lost_h2osoi[1], mc);
		  xstmp1[j][mc] = thisxs[mc];
		  if (lost_h2osoi[mc] != 0. && soil.lost0 == 0)
		    {
		      printf("\nLost 11.5 @ tracer %i t %ld\n", mc, time_var.daytime);
		      soil.lost0 = 1;
		    }
		}
	      for (mc=1; mc<=wiso.nwater+1; mc++)
		{
		  h2osoi[j  ][mc] -= thisxs[mc]/(1000.*dzsoi[j  ]);
		  h2osoi[j+1][mc] += thisxs[mc]/(1000.*dzsoi[j+1]);
		}
	    }
        }
      // nsoil
      thisxs[1] = __max(h2osoi[soil.n_soil][1]-soil.theta_s[soil.n_soil], 0.) * 1000. * dzsoi[soil.n_soil];
      xstmp1[soil.n_soil][1] = thisxs[1];
      if (thisxs[1] != 0.)
	{
	  for (mc=2; mc<=wiso.nwater+1; mc++)
	    {
#ifdef DIAG
	      if (fabs(h2osoi[soil.n_soil][1]) < 1e-15 &&
		  (h2osoi[soil.n_soil][mc] != 0. || h2osoi[soil.n_soil][1] != 0.))
		printf("\nBefore Isorat 13 tracer %i @ %ld: rare %g abundant %g\n", 
		       mc, time_var.daytime, 
		       h2osoi[soil.n_soil][mc], h2osoi[soil.n_soil][1]);
#endif
	      thisxs[mc] = thisxs[1] *
		ISORAT(&h2osoi[soil.n_soil][mc], &h2osoi[soil.n_soil][1], &lost_h2osoi[mc], &lost_h2osoi[1], mc);
	      xstmp1[soil.n_soil][mc] = thisxs[mc];
	      if (lost_h2osoi[mc] != 0. && soil.lost0 == 0)
		{
		  printf("\nLost 11.6 @ tracer %i t %ld\n", mc, time_var.daytime);
		  soil.lost0 = 1;
		}
	    }
	  for (mc=1; mc<=wiso.nwater+1; mc++)
	    {
	      h2osoi[soil.n_soil][mc] -= thisxs[mc]/(1000.*dzsoi[soil.n_soil]);
	      qdrai[mc]                   += thisxs[mc];
	    }
	}

      // sub-surface drainage (accumulate over time_var.time_step/soil.moisture_dt iterations)
      // isotopes assume ratio of deepest soil

      hkdrai    = (hk[soil.n_soil]+dhkdw[soil.n_soil]*dwat[soil.n_soil][1])*soil.moisture_dt;
      qdrai[1] += hkdrai;
      if (hkdrai != 0.)
	{
	  for (mc=2; mc<=wiso.nwater+1; mc++)
	    {
#ifdef DIAG
	      if (fabs(h2osoi[soil.n_soil][1]) < 1e-15 &&
		  (h2osoi[soil.n_soil][mc] != 0. || h2osoi[soil.n_soil][1] != 0.))
		printf("\nBefore Isorat 14 tracer %i @ %ld: rare %g abundant %g\n", 
		       mc, time_var.daytime, 
		       h2osoi[soil.n_soil][mc], h2osoi[soil.n_soil][1]);
#endif
	      rdrai = ISORAT(&h2osoi[soil.n_soil][mc], &h2osoi[soil.n_soil][1], &lost_h2osoi[mc], &lost_h2osoi[1], mc);
	      if (lost_h2osoi[mc] != 0. && soil.lost0 == 0)
		{
		  printf("\nLost 11.7 @ tracer %i t %ld\n", mc, time_var.daytime);
		  soil.lost0 = 1;
		}
	      qdrai[mc] += hkdrai*rdrai;
	    }
	}

    } // end time step iteration (iter)

#ifdef DIAG
  // Quick check for numerics
  for (j=1; j<=soil.n_soil; j++)
    {
      if (h2osoi[j][1] < 0.)
	printf("\nSoil water in layer %i < 0. after iteration: %e.\n", j, h2osoi[j][1]);
    }
#endif

  //* sub-surface drainage over time step = time_var.time_step

  for (mc=1; mc<=wiso.nwater+1; mc++)
    qdrai[mc] /= time_var.time_step;

  // limit h2osoi >= watmin. get water needed to bring h2osoi = watmin
  // from lower layer. do for all points so inner loop vectorizes
  // This is a "recharge" hack. xs is now a deficit if you will, and > 0.
  // For isotopes assume this has the isotopic ratio of the lower layer.
  
  // 1 to n_soil-1
  for (j=1; j<=soil.n_soil-1; j++)
    {
      xs[1] = __max(soil.watmin[j]-h2osoi[j][1], 0.)*1000.*dzsoi[j];
      xstmp2[j][1] = xs[1];
      if (xs[1] != 0.)
	{
	  for (mc=2; mc<=wiso.nwater+1; mc++)
	    {
#ifdef DIAG
	      if (fabs(h2osoi[j+1][1]) < 1e-15 &&
		  (h2osoi[j+1][mc] != 0. || h2osoi[j+1][1] != 0.))
		printf("\nBefore Isorat 15 tracer %i @ %ld: rare %g abundant %g\n", 
		       mc, time_var.daytime, 
		       h2osoi[j+1][mc], h2osoi[j+1][1]);
#endif
	      xs[mc] = xs[1]
		* ISORAT(&h2osoi[j+1][mc], &h2osoi[j+1][1], &lost_h2osoi[mc], &lost_h2osoi[1], mc);
	      xstmp2[j][mc] = xs[mc];
	      if (lost_h2osoi[mc] != 0. && soil.lost0 == 0)
		{
		  printf("\nLost 11.9 @ tracer %i t %ld\n", mc, time_var.daytime);
		  soil.lost0 = 1;
		}
	    }
	  for (mc=1; mc<=wiso.nwater+1; mc++)
	    {
	      h2osoi[j][mc]   += xs[mc]/(1000.*dzsoi[j]);
	      h2osoi[j+1][mc] -= xs[mc]/(1000.*dzsoi[j+1]);
	    }
	}
    }

  // n_soil
  xs[1] = __max(soil.watmin[soil.n_soil]-h2osoi[soil.n_soil][1], 0.)*1000.*dzsoi[soil.n_soil];
  xstmp2[soil.n_soil][1] = xs[1];
  if (xs[1] != 0.)
    {
      for (mc=2; mc<=wiso.nwater+1; mc++)
	{
#ifdef DIAG
	  if (fabs(h2osoi[soil.n_soil][1]) < 1e-15 &&
	      (h2osoi[soil.n_soil][mc] != 0. || h2osoi[soil.n_soil][1] != 0.))
	    printf("\nBefore Isorat 16 tracer %i @ %ld: rare %g abundant %g\n", 
		   mc, time_var.daytime, 
		   h2osoi[soil.n_soil][mc], h2osoi[soil.n_soil][1]);
#endif
	  xs[mc] = xs[1]
	    * ISORAT(&h2osoi[soil.n_soil][mc], &h2osoi[soil.n_soil][1], &lost_h2osoi[mc], &lost_h2osoi[1], mc);
	  xstmp2[soil.n_soil][mc] = xs[mc];
	  if (lost_h2osoi[mc] != 0. && soil.lost0 == 0)
	    {
	      printf("\nLost 11.10 @ tracer %i t %ld\n", mc, time_var.daytime);
	      soil.lost0 = 1;
	    }
	}
      for (mc=1; mc<=wiso.nwater+1; mc++)
	{
	  h2osoi[soil.n_soil][mc] += xs[mc]/(1000.*dzsoi[soil.n_soil]);
	  qdrai[mc]                   -= xs[mc]/time_var.time_step;
	  if (qdrai[mc] < 0.)
	    if (soil.drain0 == 0)
	      {
#ifdef DIAG
		printf("\nDrainage <0 (%e) started at time step %ld with tracer %i.\n",
		       qdrai[mc], time_var.daytime, mc);
#endif
		soil.drain0 = 1;
	      }
	}
    }

  // Check
  for (j=1; j<=soil.n_soil; j++)
    {
      if ((soil.watmin[j]-h2osoi[j][1]) > isotolerance)
	{
	  printf("\nSoil water in layer %i < watmin (%e) at time step %ld: %e (diff=%e).\n", j,
		 soil.watmin[j], time_var.daytime, h2osoi[j][1], 
		 soil.watmin[j]-h2osoi[j][1]);
	  h2osoi[j][1] = soil.watmin[j]; // shouldn't come to this
	  for (mc=2; mc<=wiso.nwater+1; mc++)
	    {
#ifdef DIAG
	      if (fabs(h2osoi[j][1]) < 1e-15 &&
		  (h2osoi[j][mc] != 0. || h2osoi[j][1] != 0.))
		printf("\nBefore Isorat 17 tracer %i @ %ld: rare %g abundant %g\n", 
		       mc, time_var.daytime, 
		       h2osoi[j][mc], h2osoi[j][1]);
#endif
	      rxs = ISORAT(&h2osoi[j][mc], &h2osoi[j][1], &lost_h2osoi[mc], &lost_h2osoi[1], mc);
	      if (lost_h2osoi[mc] != 0. && soil.lost0 == 0)
		{
		  printf("\nLost 11.11 @ tracer %i t %ld\n", mc, time_var.daytime);
		  soil.lost0 = 1;
		}
	      h2osoi[j][mc] = rxs*soil.watmin[j];
	    }
	}
    }

  // Calc total soil water and plant available water

  sm = 0.;
  sm50 = 0.;
  soil.soil_mm_root = 0.;
  for (j=1; j<=soil.n_soil; j++)
    {
      sm += h2osoi[j][1] * 1000. * dzsoi[j]; // gravel included in dzsoi
      soil.soil_mm_root += soil.root[j]*h2osoi[j][1]*1000.*dzsoi[j]; // dito
      if ((((float)((int)(10.*soil.z_soil[j])))/10.) <= 0.5) // include up to 60 cms
	sm50 += h2osoi[j][1] * 1000. * dzsoi[j];
      for (mc=1; mc<=wiso.nwater+1; mc++)
        soil.theta[j][mc] = h2osoi[j][mc];
    }
  soil.soil_mm = sm;
  soil.soil_mm_50 = sm50;

  for (mc=1; mc<=wiso.nwater+1; mc++)
    {
      soil.qdrai[mc] = qdrai[mc];
      wiso.lost[mc] = lost_h2osoi[mc];
      if (wiso.lost[mc] != 0. && soil.lost0 == 0)
	{
	  printf("\nLost 11 @ tracer %i t %ld\n",mc,time_var.daytime);
	  soil.lost0 = 1;
	}
    }

  return (!ok);
}


// ----------------------------------------------------
double ALPHA_EQU_H2O(double temp, int species)
{
  // Calculates fractionation factors during liquid-water vapour
  // equilibration at temperature temp [K]
  // species: 2: H218O; 3: HDO; else: no fractionation, return 1
  // Reference: Majoube (1971)
  double a, b, c;
  double res;

  switch (species)
    {
    case 2:
      a = +1.137e+3;
      b = -4.156e-1;
      c = -2.067e-3;
      break;
    case 3:
      a = +2.4844e+4;
      b = -7.6248e+1;
      c = +5.261e-2;
      break;
    default:
      a = 0.;
      b = 0.;
      c = 0.;
      break;
    }
  res = 1./exp((a/temp + b)/temp +c);
  return res;
}


// ----------------------------------------------------
double ALPHA_KIN_H2O(double rs, double rb, int species, int merlivat)
{
  // Calculates kinetic fractionation factor for water vapour diffusion
  // through a pure molecular diffusion resistance (rs)
  // and a boundary layer resistance (rb).
  // Assume 2/3 power law for boundary layer.
  // species: 2: H218O; 3: HDO; else: no fractionation, return 1
  // Reference: Farquhar & Lloyd (1993) - weighting
  //            Merlivat (1978) - molecular diffusision fractionation factors
  //            Cappa et al. (2003) - dito
  double alpha_k=1., res=0.;
  //
  if (merlivat == 1)
    {
      switch (species)
        {
        case 2:
          alpha_k = 0.9727;
          break;
        case 3:
          alpha_k = 0.9755;
          break;
        default:
          alpha_k = 1.;
          break;
        }
    }
  else
    {
      switch (species)
        {
        case 2:
          alpha_k = 0.9691;
          break;
        case 3:
          alpha_k = 0.9839;
          break;
        default:
          alpha_k = 1.;
          break;
        }
    }
  //
  if ((rs+rb) != 0.)
    res = (alpha_k*rs+pow(alpha_k, 2./3.)*rb) / (rs+rb);
  else
    res = 0.;
  //
  return res;
}


// ----------------------------------------------------
double DELTA(double rare, double abundant, double standard)
{
  // Calculates delta value for isotopes expressed in per mil
  double res;

  if (abundant == 0. || standard == 0.)
    res = -1.;
  else
    res = rare/abundant/standard-1.;

  return res;
}


// ----------------------------------------------------
double DELTA1000(double rare, double abundant, double standard)
{
  // Calculates delta value for isotopes expressed in per mil
  double res;

  if (abundant == 0. || standard == 0.)
    res = -1.;
  else
    res = rare/abundant/standard-1.;

  return res*1000.;
}


// ----------------------------------------------------
double DELTA_H2O(double rare, double abundant, int species)
{
  // Calculates delta value for water isotopes
  double res;

  if (abundant == 0.)
    res = -1.;
  else
    res = rare/abundant/wiso.vsmow[species]-1.;

  return res;
}


// ----------------------------------------------------
double DELTA1000_H2O(double rare, double abundant, int species)
{
  // Calculates delta value for water isotopes expressed in per mil
  double res;

  if (abundant == 0.)
    res = -1.;
  else
    res = rare/abundant/wiso.vsmow[species]-1.;

  return res*1000.;
}


// ----------------------------------------------------
double INVDELTA1000_H2O(double iso, int species)
{
  // Calculates ratio from water isotope value expressed in per mil
  return (iso/1000.+1.)*wiso.vsmow[species];
}


// ----------------------------------------------------
double INVDELTA1000(double iso)
{
  // Calculates ratio of isotope ratio to standard from delta value expressed in per mil
  return iso/1000.+1.;
}


// ----------------------------------------------------
double VTF(double T, int species)
{
  // Vogel-Tamman-Fulcher relationship: D(T)=D0*exp(-B0/(T-T0))
  // T [K]
  double A0 = 100e-09;
  double A1 = 577.;
  double A2 = 145.;
  double c1, res;
  
  switch (species)
    {
    case 2:
      c1 = sqrt(18./20.*(20.+18.)/(18.+18.));
      break;
    case 3:
      c1 = sqrt(18./19.*(19.+18.)/(18.+18.));
      break;
    default:
      c1 = 1.;
      break;
    }

  res = c1*A0*exp(-A1/(T-A2));
    
  return res;
}


// ----------------------------------------------------
int iswhite(char c)
{
  int is_white_char;
  is_white_char = (c <= 0x20 || c >= 0x7F);
  return is_white_char;
}


// ----------------------------------------------------
int find_last_nonwhite_character(char *str)
{
  int len, rpos;
  len = strlen(str);
  for (rpos=len-1; rpos >= 0; rpos--)
    if (iswhite(str[rpos]) == 0)
      break;  
  return rpos;
}


// ----------------------------------------------------
char *ltrim(char *s)
{
  while (*s == ' ' || *s == '\t')
    ++s;
  return s;
}


// ----------------------------------------------------
int LE_WISO()
{
  // Calculates latent heat of transpiration in accordance with leaf water isotopes
  //   cf. Cuntz et al. (2007)
  int j;

  int ok=1;

  // isotope independent variables
  for (j=1; j<=jtot; j++)
    { 
      // wi
      prof.sun_wi[j] = ES(prof.sun_tleaf[j]+TN0)*100./met.press_Pa; // [Pa]
      prof.shd_wi[j] = ES(prof.shd_tleaf[j]+TN0)*100./met.press_Pa; // [Pa]
      // wa
      prof.wa[j] = prof.rhov_air_filter[j][1]*(prof.tair_filter[j]+TN0)*Rw/met.press_Pa; // [Pa]
      // h
      if (prof.wa[j] < 0.)
	printf("Abs. humidity wa < 0. (%e) for layer %i at time %ld.\n", prof.wa[j], j, time_var.daytime);
      if (set_switch.no_neg_water_flux == 1)
	{
	  prof.sun_h[j] = __max(__min(prof.wa[j]/prof.sun_wi[j], 1.), 0.);
	  prof.shd_h[j] = __max(__min(prof.wa[j]/prof.shd_wi[j], 1.), 0.);
	}
      else
	{
	  prof.sun_h[j] = prof.wa[j]/prof.sun_wi[j];
	  prof.shd_h[j] = prof.wa[j]/prof.shd_wi[j];
	}

      // [m/s] -> [mol(H2O)/m2s]
      prof.rs_fact[j] = n_stomata_sides * met.air_density * 0.622; // m->mol
      // gross transpiration flux [mol(H2O)/m2s]
      prof.sun_gross[j] = prof.rs_fact[j]/(prof.sun_rs_save[j]+prof.sun_rbv[j])
	* prof.sun_wi[j];
      prof.shd_gross[j] = prof.rs_fact[j]/(prof.shd_rs_save[j]+prof.shd_rbv[j])
	* prof.shd_wi[j];

      // E [mol(H2O)/m2s]
      prof.sun_LEstoma_new[j] = prof.sun_gross[j]*(1.-prof.sun_h[j]);
      prof.shd_LEstoma_new[j] = prof.shd_gross[j]*(1.-prof.shd_h[j]);
    }

  return (!ok);
}


// ----------------------------------------------------
int LEAF_WISO()
{
  // Calculates leaf water isotopic composition after Farquhar & Cernusak (2005)
  //   with fixed leaf wtaer volume. It uses a Dongmann-style solution
  //   (Dongmann et al. 1974) for enrichment at the evaporating sites and does the
  //   'brave assumption' for bulk mesophyll water.
  // Cf. Cuntz et al. (2007)
  int mc, j;
  double Vm=18., Leff=0.020, waterc=55555.55555;

  int ok=1;

  // xylem
  for (mc=1; mc<=wiso.nwater+1; mc++)
    {
      soil.xylem[mc] = 0;
      for (j=1; j<=soil.n_soil; j++)
	soil.xylem[mc] += soil.theta[j][mc]*soil.root[j];
    }
  for (mc=2; mc<=wiso.nwater+1; mc++)
    {
#ifdef DIAG
      if (fabs(soil.xylem[1]) < 1e-15 &&
	  (soil.xylem[mc] != 0. || soil.xylem[1] != 0.))
	printf("\nBefore Isorat 17 tracer %i @ %ld: rare %g abundant %g\n", 
	       mc, time_var.daytime, 
	       soil.xylem[mc], soil.xylem[1]);
#endif
      soil.rxylem[mc] = ISORAT(&soil.xylem[mc], &soil.xylem[1],
				      &wiso.lost[mc], &wiso.lost[1], mc);
      if (wiso.lost[mc] != 0. && soil.lost0 == 0)
	{
	  printf("\nLost 02 @ tracer %i t %ld\n",mc,time_var.daytime);
	  soil.lost0 = 1;
	}
    }

  // Variables - isotope dependent
  for (mc=2; mc<=wiso.nwater+1; mc++)
    for (j=1; j<=jtot; j++)
      {
	// kinetic fractionation
	prof.sun_alpha_k[j][mc] = ALPHA_KIN_H2O(prof.sun_rs_save[j], prof.sun_rbv[j],
						       mc, wiso.merlivat);
	prof.shd_alpha_k[j][mc] = ALPHA_KIN_H2O(prof.shd_rs_save[j], prof.shd_rbv[j],
						       mc, wiso.merlivat);
	// equilibrium fractionation
	prof.sun_alpha_equ[j][mc] = ALPHA_EQU_H2O(prof.sun_tleaf[j]+TN0, mc);
	prof.shd_alpha_equ[j][mc] = ALPHA_EQU_H2O(prof.shd_tleaf[j]+TN0, mc);
	// Peclet, P
	prof.sun_peclet[j][mc] = prof.sun_LEstoma_new[j] * Leff / waterc
	  / VTF(prof.sun_tleaf[j]+TN0, mc);
	prof.shd_peclet[j][mc] = prof.shd_LEstoma_new[j] * Leff / waterc
	  / VTF(prof.shd_tleaf[j]+TN0, mc);
	// (1-exp(-P))/P
	prof.sun_fem[j][mc] = (1.-exp(-prof.sun_peclet[j][mc]))/prof.sun_peclet[j][mc];
	prof.shd_fem[j][mc] = (1.-exp(-prof.shd_peclet[j][mc]))/prof.shd_peclet[j][mc];
      }

  // Craig-Gordon
  for (mc=2; mc<=wiso.nwater+1; mc++)
    for (j=1; j<=jtot; j++)
      {
	prof.sun_craig[j][mc] = ((1.-prof.sun_h[j])*soil.rxylem[mc]/prof.sun_alpha_k[j][mc]
					+ prof.sun_h[j]*prof.rvapour[j][mc])
	  / prof.sun_alpha_equ[j][mc];
	prof.shd_craig[j][mc] = ((1.-prof.shd_h[j])*soil.rxylem[mc]/prof.shd_alpha_k[j][mc]
					+ prof.shd_h[j]*prof.rvapour[j][mc])
	  / prof.shd_alpha_equ[j][mc];
      }
  
  // Dongmann
  for (mc=2; mc<=wiso.nwater+1; mc++)
    for (j=1; j<=jtot; j++)
      {
	// evaporative site
	if (prof.sun_LEstoma_new[j] != 0.)
	  prof.sun_leafwater_e[j][mc] = prof.sun_craig[j][mc]
	    - (prof.sun_craig[j][mc] - prof.sun_leafwater_e_old[j][mc])
	    * exp(-prof.sun_alpha_k[j][mc]*prof.sun_alpha_equ[j][mc]
		  *prof.sun_gross[j]*time_var.time_step/Vm);
	else
	  prof.sun_leafwater_e[j][mc] = prof.sun_leafwater_e_old[j][mc];
	if (prof.shd_LEstoma_new[j] != 0.)
	  prof.shd_leafwater_e[j][mc] = prof.shd_craig[j][mc]
	    - (prof.shd_craig[j][mc] - prof.shd_leafwater_e_old[j][mc])
	    * exp(-prof.shd_alpha_k[j][mc]*prof.shd_alpha_equ[j][mc]
		  *prof.shd_gross[j]*time_var.time_step/Vm);
	else
	  prof.shd_leafwater_e[j][mc] = prof.shd_leafwater_e_old[j][mc];

	if (wiso.nofracleaf == 1) // nofrac = steady-state
	  {
	    prof.sun_leafwater_e[j][mc] = prof.sun_craig[j][mc];
	    prof.shd_leafwater_e[j][mc] = prof.shd_craig[j][mc];
	  }

	// bulk mesophyll
	prof.sun_leafwater[j][mc] = prof.sun_fem[j][mc] * prof.sun_leafwater_e[j][mc]
	  + (1.-prof.sun_fem[j][mc]) * soil.rxylem[mc];
	prof.shd_leafwater[j][mc] = prof.shd_fem[j][mc] * prof.shd_leafwater_e[j][mc]
	  + (1.-prof.shd_fem[j][mc]) * soil.rxylem[mc];
	// E*R_E
	if (prof.sun_LEstoma_new[j] != 0.)
	  prof.sun_trans_rtrans[j][mc] = prof.sun_alpha_k[j][mc] * prof.sun_gross[j]
	    * (prof.sun_alpha_equ[j][mc]*prof.sun_leafwater_e[j][mc]
	       - prof.sun_h[j]*prof.rvapour[j][mc]);
	else
	  prof.sun_trans_rtrans[j][mc] = 0.;
	if (prof.shd_LEstoma_new[j] != 0.)
	  prof.shd_trans_rtrans[j][mc] = prof.shd_alpha_k[j][mc] * prof.shd_gross[j]
	    * (prof.shd_alpha_equ[j][mc]*prof.shd_leafwater_e[j][mc]
	       - prof.shd_h[j]*prof.rvapour[j][mc]);
	else
	  prof.shd_trans_rtrans[j][mc] = 0.;
	// R_E
	if (prof.sun_LEstoma_new[j] != 0.)
	  prof.sun_rtrans[j][mc] = prof.sun_trans_rtrans[j][mc]/prof.sun_LEstoma_new[j];
	else
	  prof.sun_rtrans[j][mc] = 0.;
	if (prof.shd_LEstoma_new[j] != 0.)
	  prof.shd_rtrans[j][mc] = prof.shd_trans_rtrans[j][mc]/prof.shd_LEstoma_new[j];
	else
	  prof.shd_rtrans[j][mc] = 0.;
      }

  return (!ok);
}


// ----------------------------------------------------
int SOIL_FLUX_WISO()
{
  // Calculates leaf water isotopic composition after Farquhar & Cernusak (2005)
  //   with fixed leaf wtaer volume. It uses a Dongmann-style solution
  //   (Dongmann et al. 1974) for enrichment at the evaporating sites and does the
  //   'brave assumption' for bulk mesophyll water.
  // Cf. Cuntz et al. (2007)
  // water isotopes
  int mc;
  double litterevap[watersze], soilevap[watersze];
  double kv_soil=0., kv_litter=0.;
  double alpha_equ, rl, alpha_l, rs, alpha_s, ra, tmp=0.;

  int ok=1;

  // water
  soilevap[1] = flux.soilevap[1] * soil.latent;
  litterevap[1] = flux.litterevap[1] * soil.latent;

  // kvsoil is soil surface conductance to water vapor transfer (s m-1)
  kv_soil = 1. / (soil.rb+soil.rs);
  // kvsoil is litter surface conductance to water vapor transfer (s m-1)
  kv_litter = 1. / (soil.rb+soil.rl);

  for (mc=2; mc<=wiso.nwater+1; mc++)
    {
      // litter & soil
      alpha_equ = ALPHA_EQU_H2O(soil.tsfc+TN0, mc);
      // litter variables
#ifdef DIAG
      if (fabs(soil.theta_l[1]) < 1e-15 &&
	  (soil.theta_l[mc] != 0. || soil.theta_l[1] != 0.))
	printf("\nBefore Isorat 18 tracer %i @ %ld: rare %g abundant %g\n", 
	       mc, time_var.daytime, 
	       soil.theta_l[mc], soil.theta_l[1]);
#endif
      rl = ISORAT(&soil.theta_l[mc], &soil.theta_l[1],
		  &wiso.lost[mc], &wiso.lost[1], mc);
      if (wiso.lost[mc] != 0. && soil.lost0 == 0)
	{
	  printf("\nLost 04 @ tracer %i t %ld\n", mc, time_var.daytime);
	  soil.lost0 = 1;
	}
      alpha_l = ALPHA_KIN_H2O(soil.rl, soil.rb, mc, wiso.merlivat);
      // soil variables
#ifdef DIAG
      if (fabs(soil.theta[1][1]) < 1e-15 &&
	  (soil.theta[1][mc] != 0. || soil.theta[1][1] != 0.))
	printf("\nBefore Isorat 19 tracer %i @ %ld: rare %g abundant %g\n", 
	       mc, time_var.daytime, 
	       soil.theta[1][mc], soil.theta[1][1]);
#endif
      rs = ISORAT(&soil.theta[1][mc], &soil.theta[1][1],
		  &wiso.lost[mc], &wiso.lost[1], mc);
      if (wiso.lost[mc] != 0. && soil.lost0 == 0)
	{
	  printf("\nLost 05 @ tracer %i t %ld\n", mc, time_var.daytime);
	  soil.lost0 = 1;
	}
      alpha_s = ALPHA_KIN_H2O(soil.rs, soil.rb, mc, wiso.merlivat);
      // atmosphere
      ra = prof.rhov_air_filter[1][mc]/prof.rhov_air_filter_save[1];
      // litter evap
      if (wiso.nofraclitter == 1)
	{
	  if (litterevap[1] >= 0.)
	    litterevap[mc] = rl * litterevap[1];
	  else
	    litterevap[mc] = ra * litterevap[1];
	}
      else
	{
	  if (litterevap[1] != 0.)
	    litterevap[mc] = soil.c_litterevap_save*soil.lecoef*alpha_l*kv_litter
	      * (soil.rh_litter*alpha_equ*rl*ES(soil.tsfc+TN0)*100.
		 - prof.rhov_air_filter[1][mc]*(soil.T_air+TN0)*Rw);
	  else
	    litterevap[mc] = 0.;
	}
      // soil evap
      if (wiso.nofracsoil == 1)
	{
	  if (soilevap[1] >= 0.)
	    soilevap[mc] = rs * soilevap[1];
	  else
	    soilevap[mc] = ra * soilevap[1];
	}
      else
	{
	  if (soilevap[1] != 0.)
	    soilevap[mc] = soil.lecoef*alpha_s*kv_soil
	      * (soil.rh_soil*alpha_equ*rs*ES(soil.tsfc+TN0)*100.
		 - prof.rhov_air_filter[1][mc]*(soil.T_air+TN0)*Rw);
	  else
	    soilevap[mc] = 0.;
	}
    } // for (mc=2; mc<=wiso.nwater+1; mc++)

  if (soil.maxlitter == 1)
    for (mc=2; mc<=wiso.nwater+1; mc++)
      {
#ifdef DIAG
	if (fabs(litterevap[1]) < 1e-15 &&
	    (litterevap[mc] != 0. || litterevap[1] != 0.))
	  printf("\nBefore Isorat 20 tracer %i @ %ld: rare %g abundant %g\n", 
		 mc, time_var.daytime, 
		 litterevap[mc], litterevap[1]);
#endif
	rl = ISORAT(&litterevap[mc], &soil.litterevap_save,
		    &wiso.lost[mc], &wiso.lost[1], mc);
	if (wiso.lost[mc] != 0. && soil.lost0 == 0)
	  {
	    printf("\nLost 06 @ tracer %i t %ld\n", mc, time_var.daytime);
	    soil.lost0 = 1;
	  }
	litterevap[mc] = soil.maxlitterevap * rl;
      }

  // Checks
  for (mc=2; mc<=wiso.nwater+1; mc++)
    {
      if (soilevap[1] != 0.) 
	tmp = (1.-(fabs(soilevap[1])-fabs(soilevap[mc])/wiso.vsmow[mc])/fabs(soilevap[1]))*100.;
      else
	tmp = 100.;
      if (tmp > 0.1)
	if (soilevap[1] == 0. && soilevap[mc] != 0.)
	  printf("\n%i s1=0 %g %g bad\n", mc, soilevap[1], soilevap[mc]);
      if (litterevap[1] == 0.) 
	tmp = 100.;
      else
	tmp = (1.-(fabs(litterevap[1])-fabs(litterevap[mc])/wiso.vsmow[mc])/fabs(litterevap[1]))*100.;
      if (tmp > 0.1)
	if (litterevap[1] == 0. && litterevap[mc] != 0.)
	  printf("\n%i l1=0 %g %g bad\n", mc, litterevap[1], litterevap[mc]);
    }
  
  // mm -> W/m2
  for (mc=2; mc<=wiso.nwater+1; mc++)
    {
      flux.soilevap[mc] = soilevap[mc] / soil.latent;
      flux.litterevap[mc] = litterevap[mc] / soil.latent;
      flux.s_evap[mc] = flux.litterevap[mc] + flux.soilevap[mc];
    }

  return (!ok);
}


// ----------------------------------------------------
int CANOPY_FLUX_WISO()
{
  // Calculates canopy water isotope transpiration and evaporation
  //   from leaf level fluxes

  int mc, j;
  double rrcws;

  int ok=1;

  // isotope canopy transpiration
  for (mc=2; mc<=wiso.nwater+1; mc++)
    for (j=1; j<=jtot; j++)
      {
	prof.sun_LEstoma[j][mc] = prof.sun_trans_rtrans[j][mc]
	  * LAMBDA(prof.tair_filter_save[j]+TN0); // m/s -> J/m2s = W/m2
	prof.shd_LEstoma[j][mc] = prof.shd_trans_rtrans[j][mc]
	  * LAMBDA(prof.tair_filter_save[j]+TN0); // m/s -> J/m2s = W/m2
      }

  // isotope canopy evaporation
  for (mc=2; mc<=wiso.nwater+1; mc++)
    for (j=1; j<=jtot; j++)
      {
#ifdef DIAG
	if (fabs(prof.cws[j][1]) < 1e-15 &&
	    (prof.cws[j][mc] != 0. || prof.cws[j][1] != 0.))
	  printf("\nBefore Isorat 21 tracer %i @ %ld: rare %g abundant %g\n", 
		 mc, time_var.daytime, 
		 prof.cws[j][mc], prof.cws[j][1]);
#endif
	rrcws = ISORAT(&prof.cws[j][mc], &prof.cws[j][1],
		       &wiso.lost[mc], &wiso.lost[1], mc);
	if (wiso.lost[mc] != 0. && soil.lost0 == 0)
	  {
	    printf("\nLost 01.1 @ tracer %i t %ld\n",mc,time_var.daytime);
	    soil.lost0 = 1;
	  }
	// sun
	if (prof.sun_LEwet[j][1] != 0.)
	  if (prof.cws[j][1] != 0.)
	    prof.sun_LEwet[j][mc] = rrcws * prof.sun_LEwet[j][1];
	  else
	    printf("\nMulto problemo 01.1 @ tracer %i t %ld\n", mc, time_var.daytime);
	else
	  prof.sun_LEwet[j][mc] = 0.;
	// shd
	if (prof.shd_LEwet[j][1] != 0.)
	  if (prof.cws[j][1] != 0.)
	    prof.shd_LEwet[j][mc] = rrcws * prof.shd_LEwet[j][1];
	  else
	    printf("\nMulto problemo 01.2 @ tracer %i t %ld\n", mc, time_var.daytime);
	else
	  prof.shd_LEwet[j][mc] = 0.;
      }

  // isotope canopy evapotranspiration
  for (mc=2; mc<=wiso.nwater+1; mc++)
    for (j=1; j<=jtot; j++)
      prof.dLEdz[j][mc] = prof.dLAIdz[j] *
	(solar.prob_beam[j] * (prof.sun_LEstoma[j][mc]+prof.sun_LEwet[j][mc])
	 + solar.prob_shd[j] * (prof.shd_LEstoma[j][mc]+prof.shd_LEwet[j][mc]));

  // total isotope canopy transpiration and evaporation
  for (mc=2; mc<=wiso.nwater+1; mc++)
    {
      flux.c_evaporation[mc] = 0.;
      flux.c_transpiration[mc] = 0.;
      for (j=1; j<=jtot; j++)
	{
	  flux.c_evaporation[mc] += prof.dLAIdz[j] *
	    (prof.sun_LEwet[j][mc]*solar.prob_beam[j]
	     + prof.shd_LEwet[j][mc]*solar.prob_shd[j])
	    / LAMBDA(prof.tair_filter_save[j]+TN0);
	  flux.c_transpiration[mc] += prof.dLAIdz[j] *
	    (prof.sun_LEstoma[j][mc]*solar.prob_beam[j]
	     + prof.shd_LEstoma[j][mc]*solar.prob_shd[j])
	    / LAMBDA(prof.tair_filter_save[j]+TN0);
	}
      // total flux
      flux.c_evapotranspiration[mc] = 
	flux.c_evaporation[mc] + flux.c_transpiration[mc];
      flux.evapotranspiration[mc] = 
	flux.c_evapotranspiration[mc] + flux.s_evap[mc];
    }

  return (!ok);
}


// ----------------------------------------------------
double corr(int n, double x[], double y[])
{ // modified from pearsn.c of numerical recipes in c, 2nd edition, 1992
  // calcs pearsons correlation coefficient
  double TINY=1.0e-20;
  int j;
  double yt, xt;
  double syy=0., sxy=0., sxx=0., ay=0., ax=0.;
  
  for (j=0; j<n; j++)
    {
      ax += x[j];
      ay += y[j];
    }
  ax /= n;
  ay /= n;
  for (j=0; j<n; j++)
    {
      xt   = x[j]-ax;
      yt   = y[j]-ay;
      sxx += xt*xt;
      syy += yt*yt;
      sxy += xt*yt;
    }
  return sxy/(sqrt(sxx*syy)+TINY);
}


// ----------------------------------------------------
double ALPHA_KIN_CO18O(double rs, double rb, int species)
{ // Calculates kinetic fractionation factor for co18o diffusion
  //   through a pure molecular diffusion resistance (rs)
  //   and a boundary layer resistance (rb).
  // Assume 2/3 power law for boundary layer.
  // species: 2: CO18O; else: no fractionation, return 1.
  // Reference: Farquhar & Lloyd (1993)
  double alpha_k=1., res=0.;
  //
  switch (species)
    {
    case 2:
      alpha_k = 0.9912;
      break;
    default:
      alpha_k = 1.;
      break;
    }
  //
  if ((rs+rb) != 0.)
    res = (alpha_k*rs+pow(alpha_k, 2./3.)*rb) / (rs+rb);
  else
    res = 0.;
  //
  return res;
}


// ----------------------------------------------------
double ALPHA_EQU_CO2_H2O(double T, int species)
{ // Calculates the equilibrium fractionation factor 
  //   of CO2 equilibrated with water
  // Brenninkmeier et al., Isotope Geoscience, 1, 181-190, 1983
  //  T in [KELVIN]
  //  species: 2: CO18O; else: no fractionation, return 1.
  double alpha_eq=1.;
  //
  if (T > TN0)
    {
      switch (species)
	{
	case 2:
	  alpha_eq = 1. + (17604./T - 17.93)/1000.;
	  break;
	default:
	  alpha_eq = 1.;
	  break;
	}
    } // else Ice, T<T0 default=1
  return alpha_eq;
}

/*
// ----------------------------------------------------
void CO18O()
{ // Calculates co18o fluxes

  alpha_d = ALPHA_KIN_CO18O(rs, rb, 2);
  alpha_s = 1. + -7.2/1000.;
  alpha_x = 1. + -8.8/1000.;
  alpha_eq = ALPHA_EQU_CO2_H2O(t, 2);

// SOIL
      soil.T_l = soil.tsfc+TN0;
      soil.T_soil[i] = T_new_soil[i+1];

      zrs = ALPHA_EQU_CO2_H2O(t, 2) * zrssw

      // soil co18o flux in [mol(co18o)/m^2/s]
      z18fsoil = zalphas*(zrs*phet)

// INVASION
	
         zta = ptsm1m-tmelt
         zsoilstrest = max(psoilst-pplmin,0.)
         // splited evenly between REAL soil and soil pore space (porosity),
         // i.e. half of the water is absorbed by soil, half of the water stays in the pore space
         // => maximum 75%/2=37.5% of WSMAX in water filled pore space
         zthetaw = max(zsoilstrest/2., 0.)
         zthetaa = max(tporosity-zthetaw, 0.)
         zbunsen = 1.739*exp(-0.039*zta+0.000236*zta*zta)
         select case (ninv)
            case (1) // Tans (1998)
               zkappa = z23
               // 0.037*EXP() is the solubility of CO2 in H2O but there is only a chance of
               // one third that the 18O Atom will be exchanged from H2O to CO2.
               zkh  = 0.037*exp(0.118*(zta-25.))
               zkho = zkh/3.
               zda  = 0.14
            case (2) // Stern et al. (2001)
               // Ds=diffusity in soil, Da=diffusity in air.
               // Yin&Yury (1996)/Stern et al.(1999): Ds/Da = theta^2/porosity^2/3
               // They called this kappa=tortuosity.
               // But Tans (1998) defined: Ds = kappa*theta*Da
               // i.e. kappa=theta/porosity^2/3 (as in Stern et al.)
               // i.e. 0<kappa<porosity^1/3, compared to 2/3 at Tans.
               //   =0.73(porosity=0.4), =0.79(porosity=0.5), at kappa=porosity/2: 0.37 or 0.4.
               zkappa = zthetaa/tporosity**z23
               zkh    = 0.037*exp((zta-25.)*32.272/ptsm1m)
               zkho   = zkh/3.
               // ////// TEST: kh=100*kh
               // zkho   = zkho*100.
               // ////// TEST: kh=100*kh
               zda    = 0.139*((zta+273.25)/tmelt)**1.75
               //              => Tans_T compared to Stern_S: Kappa_S=1/2*Kappa_T
               //              => Piston_S=1/Sqrt(2)*Piston_T=0.7*Piston_T
            case default
               if (p_pe == p_io) then
                  write(nerr,*) 'co18o: ninv not 1 or 2.'
                  call finish ('co18o', 'ninv should be checked in inico18o1')
               endif
         end select
         zd18    = zalphamax*zda
         zpiston = sqrt(zkappa*zthetaw*zthetaa*zbunsen*zkho*zd18)/100. // /100 for cm/s -> m/s
         // CO2 invasion flux in [mol(CO2)/m^2/s]
         zfinv   = (zpiston*zca*papm1ll/ptm1ll)/R
         // CO2 equilibrated with soil water
         zrsi = eps2alpha(equifrac(ptsm1m))*zriiw
         // CO18O invasion flux in [mol(co18o)/m^2/s]
         z18finv = zfinv*(zrsi-zra)

// STEM
      ztmp2d   = spread(zrs, dim=2, ncopies=nl)
      ztmp2d   = zalphas*(ztmp2d*(pautol-prdl))
      z18fstem = sum(ztmp2d*plail, dim=2)

// ASS
      zalphadstar(:) = 0.
      do il=1, nl
         // assimilation
         ztmp2d(:,il)  = zalphad*zmmol*zgm(:,il)*(zca*zra-zcc(:,il)*zrl)
         ztmp2d2(:,il) = zrl
         ztmp1d        = 1.+tepsd/1000.*(1.-zcc(:,il)/zca)
         zalphadstar   = zalphadstar+ztmp1d*plail(:,il)
         if (lcarboanhy) then
            ztmp2d(:,il)  = ztmp2d(:,il) - zalphad*zmmol*zgm(:,il) &
                 *(1.-cafac2(1:klon,jrow))*zcc(:,il)*(ztmp1d*zra-zrl)
            ztmp2d2(:,il) = cafac2(1:klon,jrow)*zrl+(1.-cafac2(1:klon,jrow))*ztmp1d*zra
         endif
         ztmp1d = zalphad*zmmol*zgm(:,il)
         z18fla(:,il) = ztmp1d*zcc(:,il)*ztmp2d2(:,il)
         zagmrl(:,il) = ztmp1d*ztmp2d2(:,il)
         zagmcc(:,il) = ztmp1d*zcc(:,il)
         // leaf respiration
         if (nlres == 0) then
            where ((pass(:,il)-prdl(:,il)) < 0.) 
               ztmp2d(:,il) = zalphamax*zrl*(pass(:,il)-prdl(:,il))
               ztmp1d = zalphamax*zmmol*zgm(:,il)
               z18fla(:,il) = ztmp1d*zcc(:,il)*zrl
               zagmrl(:,il) = ztmp1d*zrl
               zagmcc(:,il) = ztmp1d*zcc(:,il)
            end where
         endif
      end do
      where (.not.island2) ztmp2d = 0.
      z18fass = sum(ztmp2d*plail, dim=2)
      where (zlai >= 1.0e-4) // 1e-2/30.
         zalphadstar = zalphadstar/zlai
      elsewhere
         zalphadstar = 0.
      end where

}
*/

// ----------------------------------------------------
// Fin
// ----------------------------------------------------


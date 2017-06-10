#include <stdio.h>
#include <math.h>
#include <float.h>
#include <io.h>
#include <stdlib.h>
#include <time.h>

#define PI  3.14159        // Pi
#define sze 41             // number of canopy layers plus one
#define sze3 121           // number of atmospheric layers plus one


/* =================================================================
     name="DispersionMatrix_V2_H_z_L=-0.13.C"

11-26-2002

	Dennis Baldocchi
	Ecosystem Science Division
	Department of Environmental Science, Policy and Management
	151 Hilgard Hall
	University of California, Berkeley
	Berkeley, CA
	Baldocchi@nature.berkeley.edu

-------------------------------------------
         PROGRAM     DispersionMatrix_V2_oak.c

       This version is Compiled on Microsoft C++  

       This program computes the dispersion matrix, according to
       Raupach (1989).  This model is coupled later to CANOAK
       the oak photosynthesis model, to compute CO2 profiles.

       A number of particles is released at 40 levels
       in a horizontally, homogeneous canopy.

       The vertical velocity w is computed with a Markov sequence.
       The algorithm after the approach of Thomson (1987)
       is used.

              dw(t) = a(w,z) dt + b(w,z) du

       where du is a random number with random increment with
       mean of zero, a variance equal one and a Gaussian probability distribution.

       The random number is drawn from the rejection
       technique described in Spanier and Gelbard (1969), Monte
       Carlo Principles and Neutron Transport.

       Tests show that this technique gives a mean of zero
       and variance of one--an improvement over earlier
       versions of this model.  Tests also show that it is
       faster than the central tendency method.  But we need to integrate between +/-
       5 standard deviations to get best results, rather than 3, as in V1.

       A uniform number of particles is released from each layer
       to assess the dispersion matrix.  Then linked to CANOAK
       the actual source-sink strength is computed.  We must also adjust the concentrations
       of the dispersion matrix, scaled to height defined at variable SZE, to the reference height in
       the respective field study.

       This model computes the random flight field for a 1-dimensional
       canopy and continuous release.

       Dispersion is only in the vertical and the canopy is
       assumed to be horizontally homogeneous.

      
       The system studied here is analogous to a
       volume averaged point source, stacked one atop another.
       Concentrations are computed on the principle of superposition.

       Since the canopy is horizontally homogeneous we can assume
       an infinite extent in the y direction, so the concentration
       does not vary in y and the source can be expressed in terms of
       unit y.

      
       In this program the particles released for a prescribe
       time duration--TIMEMAX [s]. Tests in the past show that TIMEMAX is different
       for tall forests and short crops.

       The simulation volume extends up to 3 * h, the crop height.
       This altitude is divided in 40 layers (DZ IHIGH).

	 For atmospheric stability effects, different Dispersion matrices are produced for different z/L values    

       


***********************
	 Debug Notes


11/2002

  Using Massman and Weil parameterizations for Tl in the canopy


5/21/2001.

	Changed the random number range to plus minus five. The random number is now
	computed with a mean of zero and a variance of 1.00. Before the variance was about
	0.96.

	Working with old code on lap top. Need to bring back linear profiles of sigma w and new
	Massman algorithm of Tl. Dave Bowling was finding bizarre kinks in conc profiles with the
	non-linear sigma w profiles. Converting the code to a structure based style.



7/28/99

	Got note for Cesscati finding error in Do loop for random number that needs fixing. We find that
	we were missing a DO.  Nancy also recommends replacing RNDDIV with RAND_MAX, so the random
	number generator is not compiler specific.

	seeding srand with a clock call.

	Having problems with non-linear sigma w parameterization. Substituting the linear terms
	from movpond for test.

	Now using sigw=a exp(b z/h) function

*/




// Declare subroutines


	double SIGMA (double Z);  // vertical profile of standard deviation of w
	double DW2DZ (double Z);  // vertical profile of variance of w
	double TL (double Z);     // vertical profile of Lagrangian Time scale
	void FILEOUT ();          // output information and data
	double RANDGEN ();        // random number generator
	double MIN(double x, double y);  // minimum of two numbers
	int Time_sec();					 // time routine


// Declare structures

		struct random_numbers
		{
			double low;     // lower limit of random number generator
			double high;    // upper limit of random number generator
			double delta;   // difference between high and low
			double random;  // random number
		} random;

		struct domain_description
        {
			long nlevel;			// number of canopy layers
			long upper_boundary;    // upper boundary, number of canopy heights
			long nlev[sze];
			
			double delta_z;        // distance between layers, m
		

        }domain;


		struct parcel_movement
        {
			double w;      // vertical velocity, w
			double term1;  // first term of Thomson Eq
			double term2;  // second term of Thomson Eq
			double term3;  // third term of Thomson Eq
			double z;      // vertical position, z
			double std_w;  // std deviation w
			double var_w;  // variance of w 
			double delta_t; // time step
			double sum_w;   // sum vertical velocities for mean
			double sum_var_rand;  // sum variance of random number
			double sum_random;    // sum random numbers
	
			double wwmean[sze3];
    		double consum[sze3];
			double conc[sze3];
			double DIJ[sze3][sze];
	

			long move;     // number of movement steps
			long partadd[sze3];
			long sumn;
		

		}parcel;


		struct turbulence_statistics
		{
			double sigma_h;   // sigma w at canopy height and above
			double sigma_zo;   // sigma w at the ground
			double del_sigma;  // derivation of sigma w with height
			double laglen;        // Lagrangian length scale
			double Tlh;        // turbluent time scale
			double sigma_sur;
		}turb;

// Declare integers and floats 


long IZF; 




time_t ltime, start, finish;


// Declare Constants

const double npart = 200000.;             // Number of parcels
const double timemax = 5000.0;           // time of parcel run


const double ustar = 1.;                //  changed: old=0.405 new=1, friction velocity (m s-1) 
const double HH = 33.;                    //  changed: canopy height (m) 
const double DD = 22.;					// changed: zero plane displacement	

main()
{


	
	// declare variables

	int I, PART, IHIGH, ILEVEL, NN, IHIGH1;

	long IT;

	double timescale;
	double xmod,zmod;



// File pointers for opening files


        FILE *fptr1;

        //   Thomson dispersion matrix  

            fptr1=fopen("c:\\local\\Canisotope\\DIJ2000O_z_L=-0.13._C","w"); 

        
		random.low=-5.;
		random.high=5.;
		random.delta=random.high-random.low;

		domain.nlevel=40; 



		// identify variables 

		domain.delta_z = HH / domain.nlevel;     // distance between vertical layers 

		domain.upper_boundary = 3;           // upper bound,three times canopy height, 99 m 
		IHIGH = domain.upper_boundary * domain.nlevel;    // number of vertical layers 

		turb.sigma_h = 1.21 * ustar;      //  sigma w > HH  */


		turb.sigma_zo = 0.0953;            //changed old=0.3 new=0.0874

		turb.sigma_sur= turb.sigma_zo * turb.sigma_h;  // sigma w at zero


		turb.del_sigma = (turb.sigma_h - turb.sigma_sur) / HH;


	



// seed random number with time 

   srand(Time_sec());

/*
***************************************************************
       Time step length (* TL) and Lagrangian Length scale
****************************************************************
*/

turb.Tlh = .3 * HH / ustar;                 /* use 0.3 instead of 0.1 as factor */

parcel.delta_t = .1 * turb.Tlh;                         /* time step */
turb.laglen = turb.sigma_h * turb.Tlh;                  /* Lagrangian length scale */


parcel.sumn = 0;

for (I = 1; I <= domain.nlevel; I++)
{
        /*
        ******************************************************
        number of particles per layer, redistribute npart #
        with height according to the relative source strength
        *****************************************************
        */

        domain.nlev[I] =(int) npart / domain.nlevel;
        parcel.sumn += domain.nlev[I];
}

/*
'       NN is the number of time increments that the particle travels
'       TIMEMAX = t [s]
'       Actual time of particle trave is TIMEMAX*TL.
'       Number of travel steps is t/dt.
*/

        NN = (int)(timemax / parcel.delta_t);

/*
'    ****************************************************************
'       Start the release of particles:
'    ****************************************************************
*/
        parcel.sum_random = 0;
        parcel.sum_w = 0.;
        parcel.sum_var_rand = 0;
        parcel.move = 0;

        IT = 1;                           /* particle counter */


        // RE1 = (double)rand()/RAND_MAX;


          /*
        assume values of a Gaussian distribution
        random numbers with a mean of zero and a variance of one.
        */

            
          IHIGH1 =IHIGH + 1;


      

                timescale = TL(HH);

/*
        Release of particles is carried out separately for each layer,
        starting with the lowest level.
*/

  for (ILEVEL = 1; ILEVEL <= domain.nlevel;ILEVEL++)
  {


/*

    ****************************************************************
       1-D case. Particles are released from each level z.  We have
       a continuous source and a horizontally homogeneous canopy.
    ****************************************************************
*/


    for(I = 1; I <=IHIGH; I++)
    parcel.consum[I] = 0;



/*
****************************************************************
      at each level NLEV(LEVEL) particles are released
****************************************************************
*/
        for (PART = 1; PART <= domain.nlev[ILEVEL]; PART++)
        {

        IT++;
        
            
	xmod=(double) IT;
	zmod=fmod(xmod,100);

	parcel.z = (double)ILEVEL * domain.delta_z;   // initial height at the level


if (zmod==0)
{
printf(" ILEVEL %6i Particle  %7i  height %f time steps %i \n",ILEVEL, IT, parcel.z, I) ;
}

// the initial vertical velocity

       

        random.random=RANDGEN();


//       vertical velocity, WW


        parcel.w = SIGMA(parcel.z) * random.random;

        /*
        number of particle movements
        */

        parcel.move += 1;

        /* compute mean and variance of w and random number */

        parcel.sum_random += random.random;
        parcel.sum_w += parcel.w;
        parcel.wwmean[ILEVEL] += parcel.w;
        parcel.sum_var_rand += random.random*random.random;



   // The particle starts its run




/*         for (I=1; I <= NN;I++) */

        IZF =(int) MIN((int)(parcel.z / domain.delta_z) + 1, IHIGH1);

        I=0;

        do
         {

       // Compute the vertical position and reflect z if it is zero
	   // Need to reflect also at top in future, but need deeper domain
	   
          parcel.z += parcel.w * parcel.delta_t;

          if(parcel.z <= 0)  // reflect particle if at ground 
          {
          parcel.w = -parcel.w;
          parcel.z = -parcel.z;
          }

          IZF = (int)MIN((int)(parcel.z / domain.delta_z) + 1, IHIGH1);

/*

    Compute the concentration of material in the controlled volume.

    Here we use the algorithm of Raupach (1989).  The ensemble average
    concentration considers the fact that we have an extensive,
    horizontally homogeneous and continuous source.  Information from
    every step from t=0 to t=T is used. It is identical to releasing a plane
    source with length x or ut, as does Wilson et al.
*/
       
        parcel.consum[IZF] += parcel.delta_t;

/*       Compute the new vertical velocity.
       Introduce the bias velocity for the case of inhomogeneous
       turbulence.  This is needed to prevent the false accumulation
       of material near the ground that would otherwise occur as
       fast air is brought downward, yet slower air at lower levels
       is less apt to leave. (see Wilson et al. Thompson etc.)
*/

                random.random= RANDGEN();

                /*
                wnew = -wold dt/Tl) + 1/2 dvarw/dz (1+w^2/varw)+
                 (2 varw dt/Tl)du
                */

                timescale = TL(parcel.z);
                parcel.std_w = SIGMA(parcel.z);
                parcel.var_w = parcel.std_w * parcel.std_w;
                parcel.term1 = -parcel.w * parcel.delta_t / timescale;
                parcel.term2 = .5 * DW2DZ(parcel.z) * (1. + (parcel.w * parcel.w) / parcel.var_w) * parcel.delta_t;
                parcel.term3 = pow((2. * parcel.var_w * parcel.delta_t / timescale),.5) * random.random;

          
                parcel.w += parcel.term1 + parcel.term2 + parcel.term3;


/*    ****************************************************************
'       STATISTICAL CHECK OF RANDOM NUMBER AND MEAN VERTICAL VELOCITY
'    ****************************************************************
*/
                /*
                number of occurences at height IZF and its
                mean vertical velocity
                */

                parcel.wwmean[IZF] += parcel.w;

                parcel.move += 1;
                parcel.sum_random += random.random;
                parcel.sum_w += parcel.w;
                parcel.sum_var_rand += random.random*random.random;

                I++;

               
                } while (I <=NN && IZF <= IHIGH); /*  NEXT I  Particle-position  and end of while */

               parcel.partadd[IZF] += 1;

        }  //  next particle 


    /*
    Introduce computation of concentration at each level.
    Use super-position principle to compute canopy
    concentration profile.
    */

        for (I = 1; I <= IHIGH; I++)
        parcel.conc[I] = parcel.consum[I];

       // Compute the dispersion matrix then reset concentrations 
		
		for (I=1; I<= IHIGH;I++)
        parcel.DIJ[I][ILEVEL] = (parcel.conc[I] - parcel.conc[IHIGH]) / (domain.delta_z * domain.nlev[ILEVEL]);
	


		for (I = 1; I<=IHIGH ;I++)
		{
		printf(" %6.2f   %6i\n", parcel.DIJ[I][ILEVEL], I);
		fprintf(fptr1,"%7.2f %6i \n", parcel.DIJ[I][ILEVEL],I);
		}


}	// Next ILEVEL


/*
'    ****************************************************************
'       Statistical check of random number: VARR=variance (should = 1),
'       SUMN = mean (should = 0), sumw# = mean vertical velocity
'       (should = 0)
'    ****************************************************************
*/

		parcel.sum_var_rand = (parcel.sum_var_rand - (pow(parcel.sum_random,2.0) / parcel.move)) / (parcel.move - 1.);

		printf("Variance Random Number %f \n",parcel.sum_var_rand);
    
		parcel.sum_w = parcel.sum_w / parcel.move;
		parcel.sum_random = (int)parcel.sum_random / parcel.move;
		fprintf(fptr1,"var random  mean random");
		fprintf(fptr1,"%7.4f  %7.4f \n",parcel.sum_var_rand,parcel.sumn);


		FILEOUT();

} // end of main 


double DW2DZ (double Z)
{

	double y, sigw2;
/*
    ****************************************************************
       COMPUTES ds2/dz FOR GIVEN s(z)
    ****************************************************************
*/

	if (Z < HH)
{

	 // y = 2.*Z*turb.del_sigma*turb.del_sigma*ustar*ustar+2.*turb.sigma_zo*turb.del_sigma*ustar; 
	 // linear model */


		// first compute derivative of sigw^2/u*^2, need to convert to ds2/dz 
       // sigw2=(0.0951/HH)*exp(4.264*Z/HH);
       // y=sigw2*ustar*ustar;

	 // first compute derivative of sigw^2/u*^2, need to convert to ds2/dz 

        sigw2=(turb.sigma_zo*turb.sigma_zo*4.9222/HH)*exp(4.9222*Z/HH);  //changed

        y=sigw2*ustar*ustar;



	}
	else
	y = 0.0;
	return y;
}


void FILEOUT()
{

printf(" \n");
printf(" %6i \n", parcel.move);
printf("mean r:  %6.2f \n ", parcel.sum_random);
printf("var. r:  %6.2f \n ", parcel.sum_var_rand);
printf("mean w:  %6.2f \n ", parcel.sum_w);

return;
}



double RANDGEN ()
{
double y, f_fnc, oper, random_1, random_2, random_sqrd;

		/*
		Produces a Random number with a mean of zero and variance of one and
		a range of data between -5 and 5.. 

		This program uses the Rejection technique.
		It can be adapted to compute random numbers for skewed and kurtotic
		distributions with the Gram Charlier distribution. 
		
		This was originally done, but it was subsequently found that using skewed
		distributions is ill posed and produces no new information. Plus our
		more recent analysis of du/dt of turbulence data, that is non Gaussian, produces
		a Gaussian distribution of du/dt. Hence, non Gaussian distributions seem not
		needed in normal practice for canopy turbulence.
		
        The Rejection Methods is based on Spanier and Gelbard (1969)
        Monte Carlo Principles and Neutron Transport Theory
		
      
        For the Gaussian Routine 
       
        Call two random numbers.  If RAND #2 is less than
		PDF=F(RAND #1) then accept RAND #1. else repeat

        note: in C rand() returns number between 0 and
        32767.  The routine from Basic relied on a random
        number between 0 and 1
        */

          do
        {
        random_1 = (double)rand()/RAND_MAX;
        random_2 = (double)rand()/RAND_MAX;


        // Value of x between high and low limits
		

        y = random.low + random_1 * random.delta;

        random_sqrd = y * y;

        /*
        PDF OF RAND

        Compute function with a Gausian distribution with mean of zero 
		and std of 1. But note the function is adjustable
		
        FFNC=EXP(-(RAND-MEAN)^2/(2*SIGRAN*SIGRAN))
        */

        oper = -random_sqrd / 2.0;
        f_fnc = exp(oper);
        }while (random_2 >= f_fnc);
        
        return y;
        }





double SIGMA (double Z)
{
double y, sigw;
/*
'    ****************************************************************
'       This function gives the s(w) value for height z
'       Use linear decrease in sigma with z/h, as Wilson et al
'       show for a corn canopy.
'    ****************************************************************
'
'   DELSIG=(SIGMAH-SIGMASUR)/HH
'
*/

if (Z < HH)
{
 
	// y = turb.sigma_zo+Z*turb.del_sigma; //linear model */


/* exponential model, computed as function of sigw/u* and z/h
   need to convert to sigma w */

     
         //sigw=turb.sigma_zo*exp(2.132 *Z/HH);
 
        
        
        sigw=turb.sigma_zo*exp(2.4611 *Z/HH);   //changed old: sigw=turb.sigma_zo*exp(2.132 *Z/HH);
        
        y=sigw*ustar;


	}
	else
	y = turb.sigma_h;

	return y;
}




double TL (double Z)
{
double y, A1;

/*
    ****************************************************************
       This function gives the variation of T(L) with height
    ****************************************************************


        adopt scaling values of Massman and Weil that adjust Tl profiles for different
        canopy strcutures

  u* Tl/h = A1 * (z-d)/h;  z > h

  u* Tl/h = A1 (1-d/h) gamma_3/ (sigmaw(z)/u*); z <= h

  A1 = A2 (sigmaw(z)/u* gamma_3)^1/2 (1-d/h)^-1/2

  gamma_3 = sigmaw(h)/u*

*/

                // factor of 2 comes from (1 - d/h)^-1/2; (1-0.75)^-1/2  
                //A1=0.6*sqrt(SIGMA(Z)/ustar*turb.sigma_h)* 2.;
            
            A1=0.6*sqrt(SIGMA(Z)/ustar*turb.sigma_h/ustar)/ sqrt(1-DD/HH);    //changed
  
        if (Z <= HH)
                {

        
        // The factor 0.25 = 1 - d/h = 1 - 0.75

                //y = A1* HH*0.25*turb.sigma_h/SIGMA(Z);
                
                y = A1* HH*(1-DD/HH)*turb.sigma_h/(SIGMA(Z)/ustar);   //changed
        }
                else
                {
                //      y = TLH * pow((Z / HH),.5);

                y = A1 * (Z-DD)/ustar;
        }
            
        return y;
        }


double MIN (double z, double x)
{
		double y;

		if(z < x)
		y = z;
		else
		y = x;
		return y;
}

int Time_sec() 
{
		time_t t;

		t = time(NULL);
		// printf("Time_sec: %d\n", t); 
		return (*gmtime(&t)).tm_sec;
}

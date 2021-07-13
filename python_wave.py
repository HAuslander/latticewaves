# INCOMPLETE -- python version of wave3d.cpp : Defines the entry point for the console application.
#



#NOTE THERE ARE TWO CONVENTIONS FOR REPRESENTING ARCS/FLOWS:
#  in(SHIFT) convention -- i.e., where the incoming amplitudes of a point (x,y,z) are shifted by one, and
#  out(SHIFT) convention -- where the outgoing amplitudes of a point are shifted by one
#  wherever possible, we will dispense with "in" and "out" until the last possible step, since it can lead
#  to too much confusion. In general, a function defined on an arc can be BACKSHIFT, ZERO, or FORWARDSHIFT,
#  so that if we're summing on the set of arcs associated with a given point (x,y,z), we can have:
#		BACKSHIFT version: sum at the points displaced by stepping BACK (one unit) along the implied direction (so that the x+ contribution will have the coordinate (x-1,y,z), etc.
#
#		ZERO version: sum over i_dir, of the arcs  associated with (x,y,z) -- i.e. going out of or pointing into
#
#		FORWARDSHIFT version:  sum at the points displaced by stepping FORWARD one unit along the implied direction
#




#include "stdafx.h"

Global_BRA_ind = 0
Global_KET_ind = 1
Global_NSQ_ind = 2


Global_UNIVERSELENGTH = 16 #MUST BE EVEN, or ELSE, PARITY IN GOING ACROSS THE AXES WILL BE BROKEN
Global_UNIVERSELENGTH_z = 1 #MUST BE EVEN (or just 1), or ELSE, PARITY IN GOING ACROSS THE AXES WILL BE BROKEN
Global_NDIM = 2 #DON'T FORGET TO CHANGE NDIM2
Global_NDIM2 = (2*NDIM)

Global_FORWARDSHIFT (+1)
Global_ZEROSHIFT (0)
Global_BACKSHIFT (-1)

#include <stdio.h>
#include <complex>

#include <math.h>
#include <string.h>

#include <vector>

Global_MAXLEN = 2500
Global_MAXIMPCORR = 50



static FILE *bigfpout;
static int MC_attempts = 0;



Global_MAXAMPLITUDE 200
#Global_FABS(x) ((x) < 0 ? (-(x)) : (x))
#Global_ABS(x) ((x) < 0 ? (-(x)) : (x))
Global_false (0)
Global_true (1)

#static int vec_1to2D[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
#static int mat_pick2[20][20];


#int **particle_descriptor	# pos1neg_n1, bra0ket1, x, y, z













Global_IM1 2147483563
Global_IM2 2147483399
Global_AM (1.0/IM1)
Global_IMM1 (IM1-1)
Global_IA1 40014
Global_IA2 40692
Global_IQ1 53668
Global_IQ2 52774
Global_IR1 12211
Global_IR2 3791
Global_NTAB 32
Global_NDIV (1+IMM1/NTAB)
Global_EPS 1.2e-7
Global_RNMX (1.0-EPS)

def ran2(seed):
{
	j = 0, k = 0

	idum2=123456789;
	iy=0;
	iv = np.array((NTAB,), dtype=int)
	temp = 0;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return (double) temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX



void nrerror2(char strr[])
{
	return;

}

















#undef EPS 

Global_NRANSI
##include "nrutil.h"
Global_ITMAX 100
Global_EPSZ 1.0e-8

Global_EPSminus 3.0e-8

Global_IM1 2147483563
Global_IM2 2147483399
Global_AM (1.0/IM1)
Global_IMM1 (IM1-1)
Global_IA1 40014
Global_IA2 40692
Global_IQ1 53668
Global_IQ2 52774
Global_IR1 12211
Global_IR2 3791
Global_NTAB 32
Global_NDIV (1+IMM1/NTAB)
Global_EPS 1.2e-7
Global_RNMX (1.0-EPS)
Global_ln2 0.693147180559945



static std::complex<double> iconst(0.0, 1.0);

static int    verbose;
static double NormConst;
static double newNormConst;
static int    nInitDeflts=0;

Global_NSTACK 50
static long istack[NSTACK + 1];
##undef NSTACK
static int iset=0;
static float gset;

class vacuumclass{ 
private:

public:
	long seed;
	long ndim;
	long universelength_x; 
	long universelength_y; 
	long universelength_z; # make this 0 for a 2-d simulation

	long nsteps;
	long t; # starts out at zero;
	int continuous_amplitudes; # if true, the amplitudes are eventually non-integer
	int exceeded_maxamp;

	# the grid will have 2*NDIM "arcs" going INTO each point (we'll try to be clear w.r.t. whether we're using INTO-point OUTFROM-point convention. 
	# 0(1) will be the positive(negative) x direction,
	# while 2(3) will be the positive y direction, and so on.

	# to access the incoming grid points, use the feeder(m_dir,i,j,k) function, so that if the outgoing arc in question is the positive-x arc 
	# grid[0][i][j][k], the feeder(0, i, j, k) function will return the amplitude in the arc in the positive x direction 
	# that started at the point (i-1,j,k). 

	#std::complex <double>  grid[2][2*NDIM][UNIVERSELENGTH][UNIVERSELENGTH][UNIVERSELENGTH]; # the dimensions are: bra/ket/n_squared, direction, x, y, z
	std::complex <double>  *****grid;
	float  ***barrier; # is either 1 or 0  or -1; if 1, the node perfectly reflects all incoming arcs; if 0, everything as if normal; if -1, the node is a black hole -- particles enter, but none leave.
		#later, you can define intermediate values, e.g. if the barrier is 0.5 at a given point, it reflects half and lets the other half go through, whereas if it's -0.5, it kills half.

	#std::complex <double>  gridpt[2][NDIM2]; # dummy grid point (independent of x,y,z) for use in updating grid



	struct particle_descriptor {
		short posneg; #+1 or -1
		short isreal; #0 or 1
		short bra0ket1; # 0 or 1
		short in_dir; #0,1...DIM2
		int	x;
		int y;
		int z;
		short empty0full1;
	} amp_chip[2*NDIM]; #


	#std::complex <double>  **gridptx;

	#double transmat[2*NDIM][2*NDIM];
	double **transmat;
	double *ranvec; 
	int *vecmix; # vecmix and ranvec will allow for us to randomly choose among the arcs at any site.
	
	double *weightvec; # this will be used in smartboson case, to allow us to preferentially
					# choose outgoing particle arcs that reduce the abssum at a given point.

	#

	#WRONG:::::due to the fact that the particles only travel along the black chess squares or the white ones (in whatever dimension) the ket grid
	#WRONG:::::will just be the bra grid displaced by one unit along the first axis (which would be x, in the case of an x,y,z 3-d grid)
	std::complex <double>  feeder(int i_bk, int i_dir, int x, int y, int z); # note this function's name is unambiguous in indicating that one steps BACK along direction implied to find amplitude
	std::complex <double>  outflow(int i_bk, int i_dir, int x, int y, int z); # note this function's name is unambiguous in indicating that one steps FORWARD along direction implied to find amplitude
	double  outflow_sumsq(int i_dir, int x, int y, int z); # note this function's name is unambiguous in indicating that one steps FORWARD along direction implied to find amplitude
	int isempty_pt(int i_bk, int i, int j, int k, int backORzeroORforwardSHIFT);
	std::complex <double>  nabs;
	std::complex <double>  nalgeb;
	std::complex <double>  brasumamp;
	std::complex <double>  ketsumamp;
	std::complex <double>  brasumabs;
	std::complex <double>  ketsumabs;
	std::complex <double>  sumampsq;

	std::complex <double>  nsqsumamp;
	std::complex <double>  nsqsumabs;
	std::complex <double>  nsqsum;


	std::vector <std::complex <double> > brasumamp_t;
	std::vector <std::complex <double> > ketsumamp_t;
	std::vector <std::complex <double> >  sumampsq_t;
	std::complex <double> sum_pt(int i_bk, int i, int j, int k, int backORzeroORforwardSHIFT); # shift param is either -1, 0, or +1
		#when invoking this function, try to add a helpful BACKSHIFT, ZEROSHIFT, or FORWARDSHIFT comment

	double realsumsq_pt(int i, int j, int k, int backORzeroORforwardSHIFT); # shift param is either -1, 0, or +1; NOTE THIS IS NOT ABSOLUTE SQUARED, but RATHER SQUARED
		#when invoking this function, try to add a helpful BACKSHIFT, ZEROSHIFT, or FORWARDSHIFT comment
	#std::vector<int> vecmix(NDIM2*NDIM2);
	#std::vector<double> ranvec(NDIM2*NDIM2);

	int negdir(int i_dir);

	void generateoutflow_from(int i_bk, int i, int j, int k, long *seed); 

	void addamp(long i_bk, long i_dir, long i, long j, long k, 
		long new_dir, std::complex <double>  amount, int backORzeroORforwardSHIFT); #backORzeroORforwardSHIFT can be -1, 0, or +1
	double realsum_pt(int i_bk, int i, int j, int k, int backORzeroORforwardSHIFT); # returns the real part of the sum
	double sumsq_pt(int i, int j, int k, int backORzeroORforwardSHIFT);
	#(double) sumsq_FORWARDSHIFT(int i, int j, int k, int backORzeroORforwardSHIFT);
	void mixvec(int *vecmix, double *ranvec, int n, long *seed);

	void measureamps();
	void printamps(FILE *fp, int whichdimensionisfixed, int leveloffixeddimension, int i_bk, int real0imag1);
	void get_shifted_coord(int incoming_dir, long x, long y, long z, long *x_in, long *y_in, long *z_in, int shift); # used in getting incoming coordinates
		#when invoking this function, try to add a helpful BACKSHIFT, or FORWARDSHIFT comment -- note that ZEROSHIFT is excluded as being trivial.

	int inplaceranked[NDIM2];
	int inplacerankedrealpos[NDIM2];
	int inplacerankedrealneg[NDIM2];
	int inplacerankedimagpos[NDIM2];
	int inplacerankedimagneg[NDIM2];
	int inplace2arcs[NDIM2*NDIM2]; 
	int inplacesigns[4]; 


	int isfermion;
	int isbosonstupid;
	vacuumclass(long in_ndim, long in_nsteps, long in_universelength_x, long in_universelength_y, long in_universelength_z);       #constructor
	~vacuumclass();    #destructor

};




vacuumclass::vacuumclass(long in_ndim, long in_nsteps, long in_universelength_x, long in_universelength_y, long in_universelength_z)
{
	long i, j, k, m, i_bk, i_dir, i1, i2, i3;

	continuous_amplitudes = false;
	isfermion = false;
	t = 0;
	ndim = in_ndim;
	nsteps = in_nsteps;
	universelength_x = in_universelength_x;
	universelength_y = in_universelength_y;
	universelength_z = in_universelength_z;
	exceeded_maxamp = true;
	/*
		double **array = new double* [sizeof(double*) * rows];
		array[0] = new double [sizeof(double) * rows * columns];
		for (int i = 1; i < rows; i++)
			array[i] = array[0] + i * columns;
	*/

	/*
    char** array_2D;
    const unsigned ROWS = 10;
    const unsigned COLUMNS = 10;
    # Allocate "main" array
    array_2D = new char*[ROWS];
    # Allocate each member of the "main" array
    for (i = 0; i < ROWS; ++i)
        array_2D[i] = new char[COLUMNS];
	*/

	#std::complex <double>  grid[2*NDIM][UNIVERSELENGTH][UNIVERSELENGTH][UNIVERSELENGTH]
	#/std::complex <double> ****grid;

	#gridptx=new std::complex <double>*[2];
    #for (i = 0; i < NDIM2; ++i)
    #    gridptx[i] = new std::complex <double>[NDIM2];

	/*
    char** array_2D;
    const unsigned ROWS = 10;
    const unsigned COLUMNS = 10;
    # Allocate "main" array
    array_2D = new char*[ROWS];
    # Allocate each member of the "main" array
    for (i = 0; i < ROWS; ++i)
        array_2D[i] = new char[COLUMNS];
	*/

	/*
    for (i = 0; i < ROWS; ++i)
        delete[] array_2D[i];

    delete[] array_2D;

    return 0;
	*/

Global_Nkinds 3 # bra, ket and nsq
    grid = new std::complex <double>****[Nkinds];

	for (i_bk = 0; i_bk < Nkinds; ++i_bk)
	{
		grid[i_bk] = new std::complex <double>***[2*ndim];
		for (i_dir = 0; i_dir < 2*ndim; ++i_dir)
		{
			grid[i_bk][i_dir] = new std::complex <double>**[universelength_x]; #HHHHHH
			for (i = 0; i < universelength_x; ++i)
			{
				grid[i_bk][i_dir][i] = new std::complex <double>*[universelength_y];
				for (j = 0; j < universelength_y; ++j)
				{
					grid[i_bk][i_dir][i][j] = new std::complex <double>[universelength_z];
					for (k = 0; k < universelength_z; ++k)
					{
						grid[i_bk][i_dir][i][j][k] = 0;
					}
				}
			}
		}
	}


    barrier = new float **[universelength_x];

	for (i = 0; i < universelength_x; ++i)
	{
		barrier[i] = new float *[universelength_y];
		for (j = 0; j < universelength_y; ++j)
		{
			barrier[i][j] = new float[universelength_z];
			for (k = 0; k < universelength_z; ++k)
			{
				barrier[i][j][k] = 0;
			}
		}
	}



    # Allocate "main" array
    transmat = new double*[2*ndim];
    # Allocate each member of the "main" array
    for (i = 0; i < 2*ndim; ++i)
        transmat[i] = new double[2*ndim];

	vecmix = new int[2*ndim];
    for (i = 0; i < 2*ndim; ++i)
		vecmix[i] = i;

	ranvec = new double[2*ndim];
    for (i = 0; i < 2*ndim; ++i)
		ranvec[i] = i;

	weightvec = new double[2*ndim];
    for (i = 0; i < 2*ndim; ++i)
		ranvec[i] = 0.0;



    # Deletion is performed in reversed order.
    # Pay special attention to the delete[] operator
    # which must be used to delete arrays (instead
    # of the "simple" delete)
    #
	/*
    for (i = 0; i < ROWS; ++i)
        delete[] array_2D[i];

    delete[] array_2D;

    return 0;
	*/

	brasumamp_t.resize(nsteps);
	ketsumamp_t.resize(nsteps);
	sumampsq_t.resize(nsteps);

	for(i=0;i<2*ndim;i++)
		for(j=0;j<2*ndim;j++)
			transmat[i][j] = 1.0/ndim;
	for(i=0;i<2*ndim;i+=2)
		transmat[i][i+1] = -(ndim-1.0)/ndim;
	for(i=1;i<2*ndim;i+=2)
		transmat[i][i-1] = -(ndim-1.0)/ndim;

	/*
	#wipe the universe clean
	for(i_dir=0; i_dir<2*ndim; i_dir++)
		for(i=0; i<UNIVERSELENGTH; i++)
			for(j=0; j<UNIVERSELENGTH; j++)
				for(k=0; k<UNIVERSELENGTH; k++)
					grid[i_dir][i][j][k]=0.0;
	*/


}

vacuumclass::~vacuumclass() #(long ndim, long nsteps, long universelength)
{
	long i, j, k, m, i_bk, i_dir, i1, i2, i3;

	for(i_bk=0; i_bk<Nkinds; i_bk++)
	{
		for(i_dir=0; i_dir<2*ndim; i_dir++)
		{
			for(i=0; i<universelength_x; i++)
			{
				for(j=0; j<universelength_y; j++)
				{
					#for(k=0; k<universelength; k++)
					#	delete[] (std::complex <double>****) grid[i_dir][i][j][k];

					delete[] grid[i_bk][i_dir][i][j];
				}
				delete[] grid[i_bk][i_dir][i];
			}
			delete[] grid[i_bk][i_dir];
		}
		delete[] grid[i_bk];
	}
	delete[] grid;

    #delete[] array_2D;
	for(i=0; i<universelength_x; i++)
	{
		for(j=0; j<universelength_y; j++)
		{
			#for(k=0; k<universelength; k++)
			#	delete[] (std::complex <double>****) grid[i_dir][i][j][k];

			delete[] barrier[i][j];
		}
		delete[] barrier[i];
	}
	delete[] barrier;




    for (i = 0; i < 2*ndim; ++i)
        delete[] transmat[i];

    delete[] transmat;
	delete[] vecmix;
	
    #for (i = 0; i < 2*ndim; ++i)
    #    delete[] gridptx[i];

	#delete[] gridptx;


}

int vacuumclass::negdir(int i_dir)
{
	if((i_dir%2)==0)
		return(i_dir+1);
	else
		return(i_dir-1);
}

std::complex <double> vacuumclass::feeder(int i_bk, int i_dir, int i, int j, int k) 
								# note i_dir is the "from" direction not the "to" direction
{	std::complex <double> temp;

	#if((i>=25 && i<=25) && (j>=25 && j<=25)  && (k>=25 && k<=25))
	#	i = i;
	#if((i>=26 && i<=26) && (j>=25 && j<=25)  && (k>=25 && k<=25))
	#	i = i;

	if(i_dir == 0)
		return (temp=(std::complex <double>) grid[i_bk][0][i==0 ? universelength_x-1 : i-1][j][k]);
	else if(i_dir == 1)
		return (temp=(std::complex <double>) grid[i_bk][1][i==universelength_x-1 ? 0 : i+1][j][k]);
	else if(i_dir == 2)
		return (temp=(std::complex <double>) grid[i_bk][2][i][j==0 ? universelength_y-1 : j-1][k]);
	else if(i_dir == 3)
		return (temp=(std::complex <double>) grid[i_bk][3][i][j==universelength_y-1 ? 0 : j+1][k]);
	else if(i_dir == 4)
		return (temp=(std::complex <double>) grid[i_bk][4][i][j][k==0 ? universelength_z-1 : k-1]);
	else if(i_dir == 5)
		return (temp=(std::complex <double>) grid[i_bk][5][i][j][k==universelength_z-1 ? 0 : k+1]);
	else
		return(0.0);
}

std::complex <double> vacuumclass::outflow(int i_bk, int i_dir, int i, int j, int k) 
								# note i_dir is the "from" direction not the "to" direction
{	std::complex <double> temp;

	#if((i>=25 && i<=25) && (j>=25 && j<=25)  && (k>=25 && k<=25))
	#	i = i;
	#if((i>=26 && i<=26) && (j>=25 && j<=25)  && (k>=25 && k<=25))
	#	i = i;

	if(i_dir == 0)
		return sum_pt(i_bk, i==universelength_x-1 ? 0 : i+1, j, k, ZEROSHIFT); #(temp=(std::complex <double>) grid[i_bk][0][i==0 ? universelength_x-1 : i-1][j][k]);
	else if(i_dir == 1)
		return sum_pt(i_bk, i==0 ? universelength_x-1 : i-1, j, k, ZEROSHIFT); # (temp=(std::complex <double>) grid[i_bk][1][i==universelength_x-1 ? 0 : i+1][j][k]);
	else if(i_dir == 2)
		return sum_pt(i_bk, i, j==universelength_y-1 ? 0 : j+1, k, ZEROSHIFT); # (temp=(std::complex <double>) grid[i_bk][2][i][j==0 ? universelength_y-1 : j-1][k]);
	else if(i_dir == 3)
		return sum_pt(i_bk, i, j==0 ? universelength_y-1 : j-1, k, ZEROSHIFT); # (temp=(std::complex <double>) grid[i_bk][3][i][j==universelength_y-1 ? 0 : j+1][k]);
	else if(i_dir == 4)
		return sum_pt(i_bk, i, j, k==universelength_z-1 ? 0 : k+1, ZEROSHIFT); # (temp=(std::complex <double>) grid[i_bk][4][i][j][k==0 ? universelength_z-1 : k-1]);
	else if(i_dir == 5)
		return sum_pt(i_bk, i, j, k==0 ? universelength_z-1 : k-1, ZEROSHIFT); # (temp=(std::complex <double>) grid[i_bk][5][i][j][k==universelength_z-1 ? 1 : k+1]);
	else
		return(0.0);
}
double vacuumclass::outflow_sumsq(int i_dir, int i, int j, int k) 
								# note i_dir is the "from" direction not the "to" direction
{	std::complex <double> temp;

	#if((i>=25 && i<=25) && (j>=25 && j<=25)  && (k>=25 && k<=25))
	#	i = i;
	#if((i>=26 && i<=26) && (j>=25 && j<=25)  && (k>=25 && k<=25))
	#	i = i;

	if(i_dir == 0)
		return sumsq_pt(i==universelength_x-1 ? 0 : i+1, j, k, ZEROSHIFT); #(temp=(std::complex <double>) grid[i_bk][0][i==0 ? universelength_x-1 : i-1][j][k]);
	else if(i_dir == 1)
		return sumsq_pt(i==0 ? universelength_x-1 : i-1, j, k, ZEROSHIFT); # (temp=(std::complex <double>) grid[i_bk][1][i==universelength_x-1 ? 0 : i+1][j][k]);
	else if(i_dir == 2)
		return sumsq_pt(i, j==universelength_y-1 ? 0 : j+1, k, ZEROSHIFT); # (temp=(std::complex <double>) grid[i_bk][2][i][j==0 ? universelength_y-1 : j-1][k]);
	else if(i_dir == 3)
		return sumsq_pt(i, j==0 ? universelength_y-1 : j-1, k, ZEROSHIFT); # (temp=(std::complex <double>) grid[i_bk][3][i][j==universelength_y-1 ? 0 : j+1][k]);
	else if(i_dir == 4)
		return sumsq_pt(i, j, k==universelength_z-1 ? 0 : k+1, ZEROSHIFT); # (temp=(std::complex <double>) grid[i_bk][4][i][j][k==0 ? universelength_z-1 : k-1]);
	else if(i_dir == 5)
		return sumsq_pt(i, j, k==0 ? universelength_z-1 : k-1, ZEROSHIFT); # (temp=(std::complex <double>) grid[i_bk][5][i][j][k==universelength_z-1 ? 1 : k+1]);
	else
		return(0.0);
}


void vacuumclass::measureamps()
{
	long i, j, k, m, i_dir, i1, i2, i3;
	std::complex <double>  ctemp;
	
	double temp_prev=0.0;
	sumampsq = 0.0;
	brasumamp =  ketsumamp  = 0.0;
	brasumabs =  ketsumabs  = 0.0;
	nsqsumabs =  nsqsumabs  = 0.0;

	for (i = 0; i < universelength_x; i++)
		for (j = 0; j < universelength_y; j++)
			for (k = 0; k < universelength_z; k++)
				for(i_dir=0;i_dir<NDIM2;i_dir++)
				{
					#if(i_dir ==0 && i == 26 && j == 25 && k == 25)
					#	i_dir = 0;
					sumampsq += grid[0][i_dir][i][j][k] * grid[1][i_dir][i][j][k]; #hhfhfhfh

					ctemp = brasumamp;
					brasumamp   += grid[0][i_dir][i][j][k]; 

					ketsumamp   += grid[1][i_dir][i][j][k]; 

					brasumabs   += fabs((double) (grid[0][i_dir][i][j][k].real())); # + iconst * fabs((double) (grid[0][i_dir][i][j][k].imag())); 

					ketsumabs   += fabs((double) (grid[1][i_dir][i][j][k].real())); # + iconst * fabs((double) (grid[1][i_dir][i][j][k].imag())); 

					nsqsumamp   += grid[NSQ_ind][i_dir][i][j][k]; 
					#if(nsqsumamp.real() != temp_prev)
					#	temp_prev = nsqsumamp.real();

					nsqsumabs   += fabs((double) (grid[2][i_dir][i][j][k].real()));


				}

}

void vacuumclass::addamp(long i_bk, long i_dir, long i, long j, long k, long new_dir, std::complex <double>  amount, int backORzeroORforwardSHIFT)
{ 
	if(backORzeroORforwardSHIFT == +1 || backORzeroORforwardSHIFT == -1) # originally the function was written with just this option
	{
		if(i_dir==0)
			grid[i_bk][new_dir][i==(universelength_x-1)?0:i+1][j][k] += amount;
		else if(i_dir==1)
			grid[i_bk][new_dir][i==0?(universelength_x-1):i-1][j][k] += amount;
		else if(i_dir==2)
			grid[i_bk][new_dir][i][j==(universelength_y-1)?0:j+1][k] += amount;
		else if(i_dir==3)
			grid[i_bk][new_dir][i][j==0?(universelength_y-1):j-1][k] += amount;
		else if(i_dir==4)
			grid[i_bk][new_dir][i][j][k==(universelength_z-1)?0:k+1] += amount;
		else if(i_dir==5)
			grid[i_bk][new_dir][i][j][k==0?(universelength_z-1):k-1] += amount;
		return;
	}
	/* # no need for this; the +1 and -1 cases are the same
	else if(backORzeroORforwardSHIFT == -1) 
	{
		if(i_dir==0)
			grid[i_bk][new_dir][i==0 ? (universelength_x-1) : i-1][j][k] += amount;
		else if(i_dir==1)
			grid[i_bk][new_dir][i==(universelength_x-1) ? 0 : i+1][j][k] += amount;
		else if(i_dir==2)
			grid[i_bk][new_dir][i][j==0 ? (universelength_y-1) : j-1][k] += amount;
		else if(i_dir==3)
			grid[i_bk][new_dir][i][j==(universelength_y-1) ? 0 : j+1][k] += amount;
		else if(i_dir==4)
			grid[i_bk][new_dir][i][j][k==0 ? (universelength_z-1) : k-1] += amount;
		else if(i_dir==5)
			grid[i_bk][new_dir][i][j][k==(universelength_z-1) ? 0 : k+1] += amount;
		return;
	}
	*/
	else # if(backORzeroORforwardSHIFT==0)
	{
		printf("Don't use addamp() with a shift param of zero -- just increment the grid itself, idiot.");
		return;
	}

}

void vacuumclass::get_shifted_coord(int incoming_dir, long x, long y, long z, long *x_in, long *y_in, long *z_in, int backORzeroORforwardSHIFT)
#again this was orig written as get_coord_BACKSHIFT(), i.e. referring to the "from" arcs 
{
	if(backORzeroORforwardSHIFT == -1 || backORzeroORforwardSHIFT == +1) # originally the function was written with just the -1 option, not the 0
	{
		if(incoming_dir==0)
		{
			*x_in = (x == 0 ? (universelength_x-1) : x-1);
			*y_in = y;
			*z_in = z;
		}
		else if(incoming_dir==2)
		{
			*x_in = x;
			*y_in = (y == 0 ? (universelength_y-1) : y-1);
			*z_in = z;
		}
		else if(incoming_dir==4)
		{
			*x_in = x;
			*y_in = y;
			*z_in = (z == 0 ? (universelength_z-1) : z-1);
		}
		else if(incoming_dir==1)
		{
			*x_in = (x == (universelength_x-1) ? 0 : x+1);
			*y_in = y;
			*z_in = z;
		}
		else if(incoming_dir==3)
		{
			*x_in = x;
			*y_in = (y == (universelength_y-1) ? 0 : y+1);
			*z_in = z;
		}
		else if(incoming_dir==5)
		{
			*x_in = x;
			*y_in = y;
			*z_in = (z == (universelength_z-1) ? 0 : z+1);
		}
		return;
	}
	/*
	else if(backORzeroORforwardSHIFT == +1) 
	{
		if(incoming_dir==0)
		{
			*x_in = (x == (universelength_x-1) ? 0 : x+1);
			*y_in = y;
			*z_in = z;
		}
		else if(incoming_dir==2)
		{
			*x_in = x;
			*y_in = (y == (universelength_y-1) ? 0 : y+1);
			*z_in = z;
		}
		else if(incoming_dir==4)
		{
			*x_in = x;
			*y_in = y;
			*z_in = (z == (universelength_z-1) ? 0 : z+1);
		}
		else if(incoming_dir==1)
		{
			*x_in = (x == 0 ? (universelength_x-1) : x-1);
			*y_in = y;
			*z_in = z;
		}
		else if(incoming_dir==3)
		{
			*x_in = x;
			*y_in = (y == 0 ? (universelength_y-1) : y-1);
			*z_in = z;
		}
		else if(incoming_dir==5)
		{
			*x_in = x;
			*y_in = y;
			*z_in = (z == 0 ? (universelength_z-1) : z-1);
		}
		return;
	}
	*/
	else # if(backORzeroORforwardSHIFT==0)
	{
		printf("Don't use get_shifted_coord() with a shift param of zero -- just increment the grid itself, idiot.");
		return;
	}
}


std::complex <double> conjugate(std::complex <double> x)
{
	std::complex <double> temp;

	temp = x.imag();
	temp *= -1;
	temp += x.real();

	return(temp);
}


double diff_sq(std::complex <double> x, std::complex <double> y)
{
	std::complex <double> temp;
	
	temp = ((x-y)*conjugate((x-y)));

	return((double) temp.real());
}



double vacuumclass::sumsq_pt(int i, int j, int k, int backORzeroORforwardSHIFT)
{
	/*
		int i_dir;
		std::complex <double> ctemp=0.0;

		for(i_dir=0;i_dir<NDIM2;i_dir++)
		{
			ctemp += grid[0][i_dir][i][j][k] * (conjugate(grid[1][i_dir][i][j][k]));
		}

		return (double) ctemp.real();
	*/
	std::complex <double> ctemp;
	int i_dir;

	if(backORzeroORforwardSHIFT == +1) # original choice of this function was +1
	{
		ctemp +=
			grid[0][1][i==(universelength_x-1)?0:i+1][j][k]*conj(grid[1][1][i==(universelength_x-1)?0:i+1][j][k]);
		ctemp +=
			grid[0][0][i==0?(universelength_x-1):i-1][j][k]*conj(grid[1][0][i==0?(universelength_x-1):i-1][j][k]);
		ctemp +=
			grid[0][3][i][j==(universelength_y-1)?0:j+1][k]*conj(grid[1][3][i][j==(universelength_y-1)?0:j+1][k]);
		ctemp +=
			grid[0][2][i][j==0?(universelength_y-1):j-1][k]*conj(grid[1][2][i][j==0?(universelength_y-1):j-1][k]);
		if(NDIM>2)
		{
			ctemp +=
				grid[0][5][i][j][k==(universelength_z-1)?0:k+1]*conj(grid[1][5][i][j][k==(universelength_z-1)?0:k+1]);
			ctemp +=
				grid[0][4][i][j][k==0?(universelength_z-1):k-1]*conj(grid[1][4][i][j][k==0?(universelength_z-1):k-1]);
		}

		return (double) ctemp.real();
	}
	
	else if(backORzeroORforwardSHIFT == -1)
	{
		ctemp +=
			grid[0][0][i==0?(universelength_x-1):i-1][j][k]*conj(grid[1][0][i==0?(universelength_x-1):i-1][j][k]);
		ctemp +=
			grid[0][1][i==(universelength_x-1)?0:i+1][j][k]*conj(grid[1][1][i==(universelength_x-1)?0:i+1][j][k]);
		ctemp +=
			grid[0][2][i][j==0?(universelength_y-1):j-1][k]*conj(grid[1][2][i][j==0?(universelength_y-1):j-1][k]);
		ctemp +=
			grid[0][3][i][j==(universelength_y-1)?0:j+1][k]*conj(grid[1][3][i][j==(universelength_y-1)?0:j+1][k]);
		if(NDIM>2)
		{
			ctemp +=
				grid[0][4][i][j][k==0?(universelength_z-1):k-1]*conj(grid[1][4][i][j][k==0?(universelength_z-1):k-1]);
			ctemp +=
				grid[0][5][i][j][k==(universelength_z-1)?0:k+1]*conj(grid[1][5][i][j][k==(universelength_z-1)?0:k+1]);
		}

		return (double) ctemp.real();
	}
	
	else # if(backORzeroORforwardSHIFT==0)
	{
		for(i_dir=0;i_dir<NDIM2;i_dir++)
			ctemp += grid[0][i_dir][i][j][k]*conj(grid[1][i_dir][i][j][k]);

		return (double) ctemp.real();
	}

}




double vacuumclass::realsumsq_pt(int i, int j, int k, int backORzeroORforwardSHIFT)
#the SHIFT indicates that we have to look one point away in order to find the outgoing point -- i.e. we're summing at all the points where the arcs are going to
{
	std::complex <double> ctemp=0.0;
	int i_dir;

	if(backORzeroORforwardSHIFT == +1 || backORzeroORforwardSHIFT == -1) # original choice of this function
	{
		ctemp +=
			grid[0][1][i==(universelength_x-1)?0:i+1][j][k]*grid[1][1][i==(universelength_x-1)?0:i+1][j][k];
		ctemp +=
			grid[0][0][i==0?(universelength_x-1):i-1][j][k]*grid[1][0][i==0?(universelength_x-1):i-1][j][k];
		ctemp +=
			grid[0][3][i][j==(universelength_y-1)?0:j+1][k]*grid[1][3][i][j==(universelength_y-1)?0:j+1][k];
		ctemp +=
			grid[0][2][i][j==0?(universelength_y-1):j-1][k]*grid[1][2][i][j==0?(universelength_y-1):j-1][k];
		if(NDIM>2)
		{
			ctemp +=
				grid[0][5][i][j][k==(universelength_z-1)?0:k+1]*grid[1][5][i][j][k==(universelength_z-1)?0:k+1];
			ctemp +=
				grid[0][4][i][j][k==0?(universelength_z-1):k-1]*grid[1][4][i][j][k==0?(universelength_z-1):k-1];
		}

		return (double) ctemp.real();
	}
	/*
	else if(backORzeroORforwardSHIFT == -1)
	{
		ctemp +=
			grid[i_bk][0][i==0?(universelength_x-1):i-1][j][k];
		ctemp +=
			grid[i_bk][1][i==(universelength_x-1)?0:i+1][j][k];
		ctemp +=
			grid[i_bk][2][i][j==0?(universelength_y-1):j-1][k];
		ctemp +=
			grid[i_bk][3][i][j==0?(universelength_y-1):j+1][k];
		if(NDIM>2)
		{
			ctemp +=
				grid[i_bk][4][i][j][k==0?(universelength_z-1):k-1];
			ctemp +=
				grid[i_bk][5][i][j][k==(universelength_z-1)?0:k+1];
		}

		return (double) ctemp.real();
	}
	*/
	else # if(backORzeroORforwardSHIFT==0)
	{
		for(i_dir=0;i_dir<NDIM2;i_dir++)
			ctemp += grid[0][i_dir][i][j][k]*grid[1][i_dir][i][j][k];
		return (double) ctemp.real();

	}

}

double vacuumclass::realsum_pt(int i_bk, int i, int j, int k, int backORzeroORforwardSHIFT)
#the SHIFT indicates that we have to look one point away in order to find the outgoing point -- i.e. we're summing at all the points where the arcs are going to
{
	std::complex <double> ctemp=0.0;
	int i_dir;

	if(backORzeroORforwardSHIFT == +1 || backORzeroORforwardSHIFT == -1) # original choice of this function
	{
		ctemp +=
			grid[i_bk][1][i==(universelength_x-1)?0:i+1][j][k];
		ctemp +=
			grid[i_bk][0][i==0?(universelength_x-1):i-1][j][k];
		ctemp +=
			grid[i_bk][3][i][j==(universelength_y-1)?0:j+1][k];
		ctemp +=
			grid[i_bk][2][i][j==0?(universelength_y-1):j-1][k];
		if(NDIM>2)
		{
			ctemp +=
				grid[i_bk][5][i][j][k==(universelength_z-1)?0:k+1];
			ctemp +=
				grid[i_bk][4][i][j][k==0?(universelength_z-1):k-1];
		}

		return (double) ctemp.real();
	}
	/*
	else if(backORzeroORforwardSHIFT == -1)
	{
		ctemp +=
			grid[i_bk][0][i==0?(universelength_x-1):i-1][j][k];
		ctemp +=
			grid[i_bk][1][i==(universelength_x-1)?0:i+1][j][k];
		ctemp +=
			grid[i_bk][2][i][j==0?(universelength_y-1):j-1][k];
		ctemp +=
			grid[i_bk][3][i][j==0?(universelength_y-1):j+1][k];
		if(NDIM>2)
		{
			ctemp +=
				grid[i_bk][4][i][j][k==0?(universelength_z-1):k-1];
			ctemp +=
				grid[i_bk][5][i][j][k==(universelength_z-1)?0:k+1];
		}

		return (double) ctemp.real();
	}
	*/
	else # if(backORzeroORforwardSHIFT==0)
	{
		for(i_dir=0;i_dir<NDIM2;i_dir++)
			ctemp += grid[i_bk][i_dir][i][j][k];
		return (double) ctemp.real();

	}

}



std::complex <double> vacuumclass::sum_pt(int i_bk, int i, int j, int k, int backORzeroORforwardSHIFT)
#the SHIFT indicates that we have to look one point away in order to find the outgoing point -- i.e. we're summing at all the points where the arcs are going to
{
	std::complex <double> ctemp=0.0;
	int i_dir;

	if(backORzeroORforwardSHIFT == +1 || backORzeroORforwardSHIFT == -1) # original choice of this function
	{
		ctemp +=
			grid[i_bk][0][i==(universelength_x-1)?0:i+1][j][k];
		ctemp +=
			grid[i_bk][1][i==0?(universelength_x-1):i-1][j][k];
		ctemp +=
			grid[i_bk][2][i][j==(universelength_y-1)?0:j+1][k];
		ctemp +=
			grid[i_bk][3][i][j==0?(universelength_y-1):j-1][k];
		if(NDIM>2)
		{
			ctemp +=
				grid[i_bk][4][i][j][k==(universelength_z-1)?0:k+1];
			ctemp +=
				grid[i_bk][5][i][j][k==0?(universelength_z-1):k-1];
		}

		return ctemp;
	}
	else if(false) #(backORzeroORforwardSHIFT == -1)
	{
		ctemp +=
			grid[i_bk][0][i==0?(universelength_x-1):i-1][j][k];
		ctemp +=
			grid[i_bk][1][i==(universelength_x-1)?0:i+1][j][k];
		ctemp +=
			grid[i_bk][2][i][j==0?(universelength_y-1):j-1][k];
		ctemp +=
			grid[i_bk][3][i][j==0?(universelength_y-1):j+1][k];
		if(NDIM>2)
		{
			ctemp +=
				grid[i_bk][4][i][j][k==0?(universelength_z-1):k-1];
			ctemp +=
				grid[i_bk][5][i][j][k==(universelength_z-1)?0:k+1];
		}

		return ctemp;
	}
	else # if(backORzeroORforwardSHIFT==0)
	{
		for(i_dir=0;i_dir<NDIM2;i_dir++)
			ctemp += grid[i_bk][i_dir][i][j][k];
		return ctemp;

	}

}

int vacuumclass::isempty_pt(int i_bk, int i, int j, int k, int backORzeroORforwardSHIFT)
#the SHIFT indicates that we have to look one point away in order to find the outgoing point -- i.e. we're summing at all the points where the arcs are going to
{
	std::complex <double> ctemp=0.0;
	int i_dir;

	if(backORzeroORforwardSHIFT == +1 ) # original choice of this function
	{
		if(grid[i_bk][0][i==(universelength_x-1)?0:i+1][j][k] != 0.0)
			return(false);
		else if(grid[i_bk][1][i==0?(universelength_x-1):i-1][j][k] != 0.0)
			return(false);
		else if(grid[i_bk][2][i][j==(universelength_y-1)?0:j+1][k] != 0.0)
			return(false);
		else if(grid[i_bk][3][i][j==0?(universelength_y-1):j-1][k] != 0.0)
			return(false);
		else if(NDIM>2)
		{
			if(grid[i_bk][4][i][j][k==(universelength_z-1)?0:k+1] != 0.0)
				return(false);
			else if(grid[i_bk][5][i][j][k==0?(universelength_z-1):k-1] != 0.0)
				return(false);
			else
				return(true);
		}

		return true;
	}
	else if(backORzeroORforwardSHIFT == -1)
	{
		if(grid[i_bk][0][i==0?(universelength_x-1):i-1][j][k] != 0.0)
			return(false);
		else if(grid[i_bk][1][i==(universelength_x-1)?0:i+1][j][k] != 0.0)
			return(false);
		else if(grid[i_bk][2][i][j==0?(universelength_y-1):j-1][k] != 0.0)
			return(false);
		else if(grid[i_bk][3][i][j==(universelength_y-1)?0:j+1][k] != 0.0) 
			return(false);
		else if(NDIM>2)
		{
			if(grid[i_bk][4][i][j][k==0?(universelength_z-1):k-1] != 0.0)
				return(false);
			else if(grid[i_bk][5][i][j][k==(universelength_z-1)?0:k+1] != 0.0)
				return(false);
			else
				return(true);
		}

		return true;
	}
	else # if(backORzeroORforwardSHIFT==0)
	{
		for(i_dir=0;i_dir<NDIM2;i_dir++)
			ctemp += grid[i_bk][i_dir][i][j][k]*conj(grid[i_bk][i_dir][i][j][k]);
		if(ctemp.real() != 0.0)
			return(false);
		else
			return(true);
	}

}


/*
std::complex <double> vacuumclass::out_sum(int i_bk, int i, int j, int k)
{
		std::complex <double> ctemp=0.0;

		ctemp +=
			grid[i_bk][0][i==(universelength_x-1)?0:i+1][j][k];
		ctemp +=
			grid[i_bk][1][i==0?(universelength_x-1):i-1][j][k];
		ctemp +=
			grid[i_bk][2][i][j==(universelength_y-1)?0:j+1][k];
		ctemp +=
			grid[i_bk][3][i][j==0?(universelength_y-1):j-1][k];
		ctemp +=
			grid[i_bk][4][i][j][k==(universelength_z-1)?0:k+1];
		ctemp +=
			grid[i_bk][5][i][j][k==0?(universelength_z-1):k-1];

		return (double) ctemp.real();

}
*/


void vacuumclass::printamps(FILE *fp, int whichdimensionisfixed, int leveloffixeddimension, int i_bk, int real0imag1)
{
	long i, j, k, m, i_dir, i1, i2, i3;
	std::complex <double>  ctemp;
	
	double temp_prev=0.0;
	sumampsq = 0.0;

	fprintf(fp,"\n\n");

	if(whichdimensionisfixed != 3)
	{	
		printf("Arg#2 in vacumclass::printamps(), 'whichdimensionisfixed', must be 3.");
		return;
	}

	for (i = 0; i < universelength_x; i++)
	{
		for (j = 0; j < universelength_y; j++)
			#for (k = 0; k < universelength_z; k++)
				#for(i_dir=0;i_dir<NDIM2;i_dir++)
				{
					#if(i==4 && j==3)
					#	i=i;
					ctemp = sum_pt(i_bk,i,j,leveloffixeddimension, ZEROSHIFT); # note we're using zeroshift
					fprintf(fp, "%f\t", (float) (real0imag1 == 0 ? real(ctemp)  : imag(ctemp)));
				}
		fprintf(fp, "\n");
	}
	fprintf(fp,"\n\n");

}














#EDIT FROM HERE 10-may-06

/*

	for (j=2;j<=n;j++) {
		a=arr[j];
		b=brr[j];
		i=j-1;
		while (i > 0 && arr[i] > a) {
			arr[i+1]=arr[i];
			brr[i+1]=brr[i];
			i--;
		}
		arr[i+1]=a;
		brr[i+1]=b;
	}
*/

void piksr2(int n, double *arr, int *brr)
{
	int i,j;
	double a;
	int b;

	for (j=1;j<n;j++) {
		a=arr[j];
		b=brr[j];
		i=j-1;
		while (i >= 0 && arr[i] > a) {
			arr[i+1]=arr[i];
			brr[i+1]=brr[i];
			i--;
		}
		arr[i+1]=a;
		brr[i+1]=b;
	}
}




void vacuumclass::mixvec(int *vecmix, double *ranvec, int n, long *seed)
{
	int i, i2, imin;
	double temp;

	for(i=0;i<n;i++)
		(ranvec[i], seed) = ran2(seed);
	piksr2(n, ranvec, vecmix);
}


void vacuumclass::generateoutflow_from(int i_bk, int i, int j, int k, long *seed)
{


	
	int	i_dir, j_dir, i_pick, i_generate, j_bk, i_sign, in_dir, out_dir;
	long idum, ltemp, i_l, i_randdir, i_randdir2, i_mod, i_rand, i_temp, i_full, ntakencareof;
	double temp, dtemp, sumwt;
	long sum_inarcs_real, sum_inarcs_imag; #, abssum_inarcs_real, abssum_inarcs_imag, mod_real, mod_imag, absmod_real, absmod_imag;
	double ran2(long *idum);
	int sumnonzero=0, isign;
	std::complex <double> ctemp = 0.0, newamp, iconstant = (0.0, 1.0);
	long x_in, y_in, z_in;

	if(barrier[i][j][k] == 1)  # then barrier reflects perfectly
	{
		for(i_dir=0;i_dir<NDIM2;i_dir++)
		{
			grid[i_bk][i_dir][i][j][k] = feeder(i_bk, negdir(i_dir), i, j, k);
		}
		return;
	}
	else if(barrier[i][j][k] == -1) # then barrier absorbs perfectly
	{
		for(i_dir=0;i_dir<NDIM2;i_dir++)
		{
			grid[i_bk][i_dir][i][j][k] = 0;
		}
		return;
	}



	if((UNIVERSELENGTH % 2 == 1) || (UNIVERSELENGTH_z % 2 == 1 && UNIVERSELENGTH_z != 1))
	{
		printf("universe lengths must be even for this to work, otherwise at periphery, you will have cells with odd parity influencing one another. \n");
		exit;
	}
		

	if(true ||  continuous_amplitudes==true)
	{
		for(i_dir=0;i_dir<NDIM2;i_dir++)
		{
			#find the sum of the incoming amplitudes
			ctemp = 0.0;
			for(j_dir=0;j_dir<NDIM2;j_dir++)
			{
				ctemp += feeder(i_bk, j_dir, i, j, k)*transmat[i_dir][j_dir];
			}
			grid[i_bk][i_dir][i][j][k] = ctemp;
		}	
	}
	else if(false && isbosonstupid == true) # by this we mean that all particles are individual particle-generators, each spanning one back-particle and two random particles
	{

		for(i_dir=0;i_dir<NDIM2;i_dir++)
		{
			if(feeder(i_bk, i_dir, i, j, k) == 0.0)
				continue;
			for(j_bk=0; j_bk<2; j_bk++)  # NOT bra and ket, but real and imag
			{
				#added PREAMBLE in MAY-07:
				#randomly select largest possible group of DIM (i.e. 3 or 2) particles in any set of incoming nodes
				#and treat the resultant sum as a continuous set of amplitudes
				#as for the leftovers that were not able to be put into sets of DIM
				#continue as before, i.e.:

				#for each incoming particle, generate an integer (probabilistic) outflow
				#via the following easy method -- generate a flow of same magnitude in the 
				#opposite direction and then pick two other DISTINCT arcs (including possible the
				#opposite direction) to put a partcle of same sign (annihilate oppositely
				#charged particles in opposite direction if necessary
				#NOTE: None of the earlier versions of this routine (this was done in Apr-07) allow for 2-dim simulations

				#test to see if you've maxed out the amplitude
				#test to see if you've maxed out the amplitude
				ltemp = (long) (j_bk==0 ? feeder(i_bk, i_dir, i, j, k).real() : feeder(i_bk, i_dir, i, j, k).imag());
				if(ltemp > MAXAMPLITUDE)
					ltemp = ltemp;

				#test to see if incoming amplitude is zero -- no need to continue if it is.
				if(ltemp==0)
					continue;
				# if we've made it here, we have to generate outflow
				for(i_l=(long)fabs((double) ltemp);i_l>0;i_l--)
				{
					if(j_bk==0)
					{
						grid[i_bk][negdir(i_dir)][i][j][k] += (-DSIGN(ltemp));

						/*
						for(i_randdir=0;i_randdir<NDIM2;i_randdir++)
						{
 							grid[i_bk][i_randdir][i][j][k] += (DSIGN(ltemp))/NDIM;
						}
						*/

						
						(i_randdir, seed) = int( (ran2(seed)*NDIM2) )
						grid[i_bk][i_randdir][i][j][k] += DSIGN(ltemp)
						(i_randdir, seed) = (int (ran2(seed)*(NDIM2))) 
						#if(i_randdir2 >= i_randdir)
						#	i_randdir2++; # now we've picked a random direction, but distinct from the first one.
						grid[i_bk][i_randdir2][i][j][k] += DSIGN(ltemp);
						
					}
					else
					{
						grid[i_bk][negdir(i_dir)][i][j][k] += int(np.sign(ltemp)*(-iconstant));
						(i_randdir, seed) = (int (ran2(seed)*(NDIM2))) 
						grid[i_bk][i_randdir][i][j][k] += DSIGN(ltemp);
						(i_randdir2, seed) = (int (ran2(seed)*(NDIM2))) 
						#if(i_randdir2 >= i_randdir)
						#	i_randdir2++; # now we've picked a random direction, but distinct from the first one.
						grid[i_bk][i_randdir2][i][j][k] += DSIGN(ltemp);
					}
				}
			}
		}


	}


	else if(isfermion == true)
	{
		# first find how many of the incoming flows have nodes
		
		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS

		#VERY IMPORTANT -- MAKE SURE THAT WHEN YOU PICK ONE NEGATIVE BACKFLOW
		#AND TWO POSITIVE ONES, THAT THE TWO POSITIVE ONES ARE IN DIFFERENT ARCS


		for(i_dir=0;i_dir<NDIM2;i_dir++)
		{
			if(feeder(i_bk, i_dir, i, j, k) != (std::complex <double>) 0.0)
			{
				inplaceranked[sumnonzero] = i_dir;
				sumnonzero++;
			}
		}
		if(sumnonzero==0)
		{
			for(i_dir=0;i_dir<NDIM2;i_dir++)
			{
				grid[i_bk][i_dir][i][j][k]=0.0;
				return;
			}
		}

		for(i_dir=0;i_dir<NDIM2;i_dir++)
		{
			#pick one of the inflows
			temp=ran2(seed);
			i_pick = inplaceranked[(int) ((temp)==1.0 ? sumnonzero-1 : temp*sumnonzero)];
			# all the flows in i_dir are attributable to amplitude in i_pick times sumnonzero
			# 
			#temp = real(grid[i_bk][i_dir][i][j][k]);
			grid[i_bk][i_dir][i][j][k]=0.0;


			#Do real, then imaginary
			#first, the real
			for(j_bk=0;j_bk<2;j_bk++) # NOT bra and ket, but real and imag
			{
				temp = (j_bk==0 ? real(feeder(i_bk, i_pick, i, j, k)) : imag(feeder(i_bk, i_pick, i, j, k))); #this is an integer
				ltemp = sumnonzero*((int) (fabs(temp)));
				if(ltemp==0)
					continue;
				if(ltemp > MAXAMPLITUDE)
				{
					exceeded_maxamp = true;
					#return;
				}
				isign = DSIGN(temp)*DSIGN(transmat[i_pick][i_dir]);
				for(i_generate=0;i_generate<ltemp;i_generate++)
				{
					dtemp = ran2(seed);
					if(j_bk==0)
						grid[i_bk][i_dir][i][j][k] += isign*((int) (dtemp>fabs(transmat[i_pick][i_dir]) ? 0.0 : 1.0));
					else
					{
						ctemp = isign*((int) (dtemp>fabs(transmat[i_pick][i_dir]) ? 0 : 1.0));
						ctemp *= iconstant;
						grid[i_bk][i_dir][i][j][k] += ctemp;

					}
				}
			}
		}
	}
	else #is boson -- if abssum of particles coming into a node is not multiple of NDIM, pick out  up to NDIM stragglers.
		#treat the rest as in the continuous amplitude case; treat the stragglers as if they were stupidbosons (see above)
	{
		#first, wipe clean...
		for(i_dir=0;i_dir<NDIM2;i_dir++)
			grid[i_bk][i_dir][i][j][k]=0.0;

		#the amp_chip array will hold the info related to the stragglers.
		for(i_dir=0;i_dir<NDIM2;i_dir++)
		{
			amp_chip[i_dir].empty0full1=0;
		}
		#added PREAMBLE in MAY-07:
		#again, we randomly select largest possible group of DIM (i.e. 3 or 2) particles in any set of incoming nodes
		#and treat the resultant sum as a continuous set of amplitudes
		#as for the leftovers that were not able to be put into sets of DIM


		mixvec(vecmix, ranvec, NDIM2, seed);
		#do real first
		sum_inarcs_real = 0;
		i_temp=0;
		i_mod=0;
		do
		{
			dtemp=0;
			for(i_dir=0;i_dir<NDIM2;i_dir++)
			{
				dtemp += (float) feeder(i_bk,i_dir,i,j,k).real();# (float) grid[i_bk][i_dir][x_in][y_in][z_in].real();
				sum_inarcs_real = (long) (dtemp + .0001*DSIGN(dtemp));
			}
			sum_inarcs_real = (long) (dtemp + .0001*DSIGN(dtemp));

			
			if(sum_inarcs_real % NDIM != 0) # if this is true, must find first nonzero particle
			{
				for(i_dir=0;i_dir<NDIM2;i_dir++)
				{
					dtemp = (float) feeder(i_bk,vecmix[i_dir],i,j,k).real();
					i_sign = DSIGN(dtemp);

					if(dtemp != 0.0 && i_sign == DSIGN(sum_inarcs_real % NDIM))
					{
						get_shifted_coord(vecmix[i_dir], i, j, k, &x_in, &y_in, &z_in, BACKSHIFT);
						#i_sign = DSIGN(grid[i_bk][vecmix[i_dir]][x_in][y_in][z_in].real());

						amp_chip[i_mod].posneg = i_sign;
						amp_chip[i_mod].isreal = 1; 
						amp_chip[i_mod].bra0ket1 = i_bk;
						amp_chip[i_mod].in_dir = vecmix[i_dir];
						amp_chip[i_mod].x = x_in;
						amp_chip[i_mod].y = y_in;
						amp_chip[i_mod].z = z_in;
						amp_chip[i_mod].empty0full1=1;
						grid[i_bk][vecmix[i_dir]][x_in][y_in][z_in] -= i_sign; #decrement the amplitude
						sum_inarcs_real -= i_sign;
						i_mod++;
						i_dir=NDIM2; #break loop
					}
				}
			}
			i_temp++;
		} while(sum_inarcs_real % NDIM != 0);

		if(i_temp>NDIM2)
		{
			printf("FATAL ERROR: screwed up modulo count...\n");
			return;
		}
		#now do imaginary part
		mixvec(vecmix, ranvec, NDIM2, seed);
		i_mod=0;
		i_temp=0;
		do
		{
			dtemp=0;
			for(i_dir=0;i_dir<NDIM2;i_dir++)
			{
				dtemp += (float) feeder(i_bk,i_dir,i,j,k).imag();# (float) grid[i_bk][i_dir][x_in][y_in][z_in].real();
			}
			sum_inarcs_imag = (long) (dtemp + .0001*DSIGN(dtemp));

			if(sum_inarcs_imag  % NDIM != 0) # if this is true, must find first nonzero particle
			{
				for(i_dir=0;i_dir<NDIM2;i_dir++)
				{
					dtemp = (float) feeder(i_bk,vecmix[i_dir],i,j,k).imag();
					i_sign = DSIGN(dtemp);

					if(dtemp != 0.0 && i_sign == DSIGN(sum_inarcs_imag % NDIM))
					{
						get_shifted_coord(vecmix[i_dir], i, j, k, &x_in, &y_in, &z_in, BACKSHIFT);
						#i_sign = DSIGN(grid[i_bk][vecmix[i_dir]][x_in][y_in][z_in].real());
						get_shifted_coord(vecmix[i_dir], i, j, k, &x_in, &y_in, &z_in, BACKSHIFT);
						i_sign = DSIGN(grid[i_bk][vecmix[i_dir]][x_in][y_in][z_in].imag());

						amp_chip[i_mod].posneg = i_sign;
						amp_chip[i_mod].isreal = 0; # note this line is different in previous pass
						amp_chip[i_mod].bra0ket1 = i_bk;
						amp_chip[i_mod].in_dir = vecmix[i_dir];
						amp_chip[i_mod].x = x_in;
						amp_chip[i_mod].y = y_in;
						amp_chip[i_mod].z = z_in;
						amp_chip[i_mod].empty0full1=1;
						ctemp = iconstant * ((double) (-1. * i_sign)); # note this line is different  in previous pass
						grid[i_bk][vecmix[i_dir]][x_in][y_in][z_in] += ctemp; #decrement the amplitude # note this line is different in previous pass
						sum_inarcs_real -= i_sign;
						i_mod++;
						i_dir=NDIM2; #break loop
					}
				}
			}
			i_temp++;
		} while(sum_inarcs_imag % NDIM != 0);

		if(i_temp>NDIM2)
		{
			printf("FATAL ERROR: screwed up modulo count...\n");
			return;
		}



		#we know that the net output in any arc will be integral (given that the amplitude sum is a multiple of DIM)
		#so that now , treat the large (multiple of NDIM) set as if it were CONTINUOUS AMPLITUDES; you can do this because you are assured the output will be integer

		# add amplitudes (a loop over all directions, then also a special contribution for the backflow)
		for(i_dir=0;i_dir<NDIM2;i_dir++)
		{
			#find the sum of the incoming amplitudes
			ctemp = 0.0;
			for(j_dir=0;j_dir<NDIM2;j_dir++)
			{
				ctemp += feeder(i_bk, j_dir, i, j, k)*transmat[i_dir][j_dir];
			}
			if(ctemp.real() != 0.0)
				j_dir = j_dir;
			#now add a little error term (of like magnitude) to both real and imaginary before truncating them;
			#if this routine did integer arithmetic, this fudging would be unnecessary, as all the outputs from this loop are complex integers.
			#grid[i_bk][i_dir][i][j][k] = ctemp;
			grid[i_bk][i_dir][i][j][k] += (long) (ctemp.real() + .0001*DSIGN(ctemp.real())); # this will keep everything integral, at least for the reals.
			ltemp = ((long) (ctemp.imag() + .0001*DSIGN(ctemp.imag())));
			grid[i_bk][i_dir][i][j][k] +=  ((double) ltemp) * iconstant;
			
		}


Global_SEA_WEIGHT_PARAM 256



		#treat the modulo bits in amp_chip as randomly generated particles -- only if amp_chip[...].empy0full1 = 1
		#for each one, generate 2 outgoing particles and a backflow particle
		#FIRST PASS, MAY-07: output arcs are randomly chosen
		for(i_pick=0;i_pick<NDIM2;i_pick++)
		{
			#OLD way of picking random direction to put in 1st (of 2) new particles (not including backflow particle)
			#i_randdir2 = ran2(seed)*NDIM2;






			#
			#
			#
			#
			#begin: NEW way of picking random direction to put in 1st (of 2) new particles (not including backflow particle)
			#
			#
			if(amp_chip[i_pick].empty0full1 == 0)
				continue;
			#1st particle:
			sumwt = 0.0;

			for(i_dir=0;i_dir<NDIM2;i_dir++)
			{
				if(amp_chip[i_pick].isreal==true)
					i_sign = DSIGN(grid[i_bk][i_dir][i][j][k].real());
				else
					i_sign = DSIGN(grid[i_bk][i_dir][i][j][k].imag());

				if(i_sign == -amp_chip[i_pick].posneg)
				{
					if(amp_chip[i_pick].isreal==true)
						weightvec[i_dir] = SEA_WEIGHT_PARAM + fabs(grid[i_bk][i_dir][i][j][k].real()); # if you replace this with weightvec[i_dir] = 1.0, you'll recover stupidboson case
					else
						weightvec[i_dir] = SEA_WEIGHT_PARAM + fabs(grid[i_bk][i_dir][i][j][k].imag()); # if you replace this with weightvec[i_dir] = 1.0, you'll recover stupidboson case
				}
				else
				{
					#weightvec[i_dir] = #SEA_WEIGHT_PARAM; # if you want more dampening, change to: weightvec[i_dir] = FMAX(0.0, SEA_WEIGHT_PARAM  - fabs(grid[i_bk][i_dir][i][j][k].real/imag()))
					if(amp_chip[i_pick].isreal==true)
						weightvec[i_dir] = FMAX(0.0, SEA_WEIGHT_PARAM - fabs(grid[i_bk][i_dir][i][j][k].real())); # if you replace this with weightvec[i_dir] = 1.0, you'll recover stupidboson case
					else
						weightvec[i_dir] = FMAX(0.0, SEA_WEIGHT_PARAM - fabs(grid[i_bk][i_dir][i][j][k].imag())); # if you replace this with weightvec[i_dir] = 1.0, you'll recover stupidboson case
				}			
					
				sumwt += weightvec[i_dir];
			}
			if(sumwt <= 0.0)
			{
				for(i_dir=0;i_dir<NDIM2;i_dir++)
					weightvec[i_dir] = 1.0/NDIM2; #normalize
			}
			else
			{
				for(i_dir=0;i_dir<NDIM2;i_dir++)
					weightvec[i_dir] /= sumwt; #normalize
			}

			dtemp = ran2(seed);
			sumwt=0.0;
			for(i_dir=0;i_dir<NDIM2;i_dir++)
			{
				sumwt += weightvec[i_dir];
				if(sumwt >= dtemp)
				{
					i_randdir = i_dir; # the original version just had: #i_randdir = ran2(seed)*NDIM2;
					i_dir = NDIM2; # break loop
				}
			}
			if(sumwt < dtemp)
				i_randdir = NDIM2-1; # take care of the case where roundoff causes dtemp > sumwt even after the loop.
			# now repeat for second particle
			dtemp = ran2(seed);
			sumwt=0.0;
			for(i_dir=0;i_dir<NDIM2;i_dir++)
			{
				sumwt += weightvec[i_dir];
				if(sumwt >= dtemp)
				{
					i_randdir2 = i_dir;  # the original version just had: #i_randdir2 = ran2(seed)*NDIM2;
					i_dir = NDIM2; # break loop
				}
			}
			if(sumwt < dtemp)
				i_randdir2 = NDIM2-1; # take care of the case where roundoff causes dtemp > sumwt even after the loop.
			#
			#
			#end: NEW way of picking random direction to put in 1st (of 2) new particles (not including backflow particle)
			#
			#
			#
			#


			#(in future implementations, try two or three directions, and preferentially choose the one resulting
			# in lowest abs magnitudes, with the odds of choosing the lower one proportional to how high the amplitude is

			grid[i_bk][i_randdir][i][j][k] += (amp_chip[i_pick].isreal==false ? iconstant : 1) * (double) amp_chip[i_pick].posneg;
			#2nd particle:
			#i_randdir2 = ran2(seed)*NDIM2;
			grid[i_bk][i_randdir2][i][j][k] += (amp_chip[i_pick].isreal==false ? iconstant : 1) * (double) amp_chip[i_pick].posneg;
			#addamp(i_bk, amp_chip[i_pick].in_dir, i, j, k, i_randdir2, (amp_chip[i_pick].isreal==false ? iconstant : 1) * (double) amp_chip[i_pick].posneg, FORWARDSHIFT);

			#backfow 
			#addamp(i_bk, amp_chip[i_pick].in_dir, i, j, k, negdir(amp_chip[i_pick].in_dir), (amp_chip[i_pick].isreal==false ? iconstant * (double) -1.0 : -1) * (double) amp_chip[i_pick].posneg, FORWARDSHIFT); # note it's a negative flow
			grid[i_bk][negdir(amp_chip[i_pick].in_dir)][i][j][k] += (amp_chip[i_pick].isreal==false ? iconstant * (double) -1.0 : -1) * (double) amp_chip[i_pick].posneg;
		}



		/*

		#old stuff (PRE-MAY-07)
		
		for(i_dir=0;i_dir<NDIM2;i_dir++)
		{
			if(grid[i_bk][i_dir][i][j][k] == 0.0)
				continue;
			for(j_bk=0; j_bk<2; j_bk++)  # NOT bra and ket, but real and imag
			{
				#added PREAMBLE in MAY-07:
				#randomly select largest possible group of DIM (i.e. 3 or 2) particles in any set of incoming nodes
				#and treat the resultant sum as a continuous set of amplitudes
				#as for the leftovers that were not able to be put into sets of DIM
				#continue as before, i.e.:

				#for each incoming particle, generate an integer (probabilistic) outflow
				#via the following easy method -- generate a flow of same magnitude in the 
				#opposite direction and then pick two other DISTINCT arcs (including possible the
				#opposite direction) to put a partcle of same sign (annihilate oppositely
				#charged particles in opposite direction if necessary
				#NOTE: None of the earlier versions of this routine (this was done in Apr-07) allow for 2-dim simulations

				#test to see if you've maxed out the amplitude
				ltemp = (long) (j_bk==0 ? grid[i_bk][i_dir][i][j][k].real() : grid[i_bk][i_dir][i][j][k].imag());
				if(ltemp > MAXAMPLITUDE)
					ltemp = ltemp;

				#test to see if incoming amplitude is zero -- no need to continue if it is.
				if(ltemp==0)
					continue;
				# if we've made it here, we have to generate outflow
				for(i_l=(long)fabs((double) ltemp);i_l>0;i_l--)
				{
					if(j_bk==0)
					{
						addamp(i_bk, i_dir, i, j, k, negdir(i_dir), -DSIGN(ltemp)); # add the mandatory backflow
						i_randdir = (int) (ran2(seed)*NDIM2);
						addamp(i_bk, i_dir, i, j, k, i_randdir, DSIGN(ltemp));
						i_randdir2 = (int) (ran2(seed)*(NDIM2)); 
						#if(i_randdir2 >= i_randdir)
						#	i_randdir2++; # now we've picked a random direction, but distinct from the first one.
						addamp(i_bk, i_dir, i, j, k, i_randdir2, DSIGN(ltemp));
					}
					else
					{
						addamp(i_bk, i_dir, i, j, k, negdir(i_dir), DSIGN(ltemp)*(-iconstant));
						i_randdir = (int) (ran2(seed)*NDIM2);
						addamp(i_bk, i_dir, i, j, k, i_randdir, DSIGN(ltemp)*iconstant);
						i_randdir2 = (int) (ran2(seed)*(NDIM2)); 
						#if(i_randdir2 >= i_randdir)
						#	i_randdir2++; # now we've picked a random direction, but distinct from the first one.
						addamp(i_bk, i_dir, i, j, k, i_randdir2, DSIGN(ltemp)*iconstant);
					}
				}
			}
		}
		*/
		#now wipe the generating particles so they'll be clean for the next step.
		#for(i_dir=0;i_dir<NDIM2;i_dir++)
		#{
		#	grid[i_bk][i_dir][i][j][k] = 0.0;
		#}


	}
}






int invec(double *vecx, long j, long ncols)
{
	long i;

	for(i=1;i<=ncols;i++)
		if(vecx[i-1] == j)
			return(1);
	return(0);
}


long getnthunused(double *vecx, long iget, long howmanycols)
{

	long i, j, k, iguess=vecx[howmanycols]; #iget-1];
	long temp=0;

	for(i=1;i<=iget+temp;i++)
	{
		if(invec(vecx, i/* +temp */ , howmanycols))
		{
			temp++;
			iguess++;
		}
	}
	return(iguess);
}


long myround(double x) # if fractional part is 0.5, round towards 0
{
	long longx = (long) x;
	if(x < 0)
	{
		if(longx - x > 0.5)
			return(longx-1);
		else
			return(longx);
	}
	else
	{
		if(x - longx > 0.5)
			return(longx+1);
		else
			return(longx);
	}
}




unsigned int gcd(unsigned int n, unsigned int m) { 
   if (n == 0) return m; 
   if (m == 0) return n; 
   while (m != n) { 
      if (n > m) n = n - m; 
      else m = m - n; 
   } 
   return n; 
} 








void getampcomponents(long n_bra, long n_ket, long n_sq, double *vecin_arr, int *isOK)
{

	/*

		Must solve the following system of equations (i.e. find a,b,c,d)

		a + b = n_bra
		c + d = n_ket
		ac - bd = n_sq

		WE WANT b AND d TO BE AS LOW IN MAGNITUDE AS POSSIBLE WHICH IMPLIES THAT
			a = ROUND(n_sq/n_k) if n_k != 0
		and
			c = ROUND(n_sq/n_b) if n_b != 0 
		now if n_b=0 then a * n_k = n_sq, so that if a is not an integer, the sol'n fails and we'll need 
		jiggle n_b by plus/minus 1 to make it work. 
		likewise, if n_k=0 then c * n_b = n_sq, so that we'll need to jiggle n_b by plus/minus 1 to make it work.
   	*/

	long absfrac, a, b, c, d, a2, b2, c2, d2, a3, b3, c3, d3, a4, b4, c4, d4;
	a = b = c = d = a2 = b2 = c2 = d2 = a3 = b3 = c3 = d3 = a4 = b4 = c4 = d4 = 0;

	if(n_bra==0)
	{
		if(n_ket==0 && n_sq == 0)
		{
			a = b = c = d = a2 = b2 = c2 = d2 = a3 = b3 = c3 = d3 = a4 = b4 = c4 = d4 = 0;
		}
		else if(n_ket==0 && n_sq != 0)
		{
			*isOK = false;
			return;
		}
		else if(n_sq % n_ket != 0)
		{
			*isOK = false;
			return;
		}
		else
		{
			a = n_sq/n_ket;
			b = -a;
			if(true) #(n_sq > 0)
			{
				c = n_ket;
				d = 0;
			/*
			}
			else
			{
			*/
				a2=a;
				b2=b;
				d2 = n_ket;
				c2 = 0;
				a3=b3=c3=d3=a4=b4=c4=d4=1;
			}
			
		}
	}
	else if(n_ket==0)
	{
		if(n_bra==0 && n_sq != 0)
		{
			*isOK = false;
			return;
		}
		else if(n_sq % n_bra != 0)
		{
			*isOK = false;
			return;
		}
		else
		{
			c = n_sq/n_bra;
			d = -c;
			if(true) #(n_sq > 0)
			{
				a = n_bra;
				b = 0;
			/*
			}
			else
			{
			*/
				b2 = n_bra;
				a2 = 0;
				c2=c;
				d2=d;
				a3=b3=c3=d3=a4=b4=c4=d4=1;
			}
		}
	}
	else
	{

		double temp;
		long gcd_bk, startpt, i_inc;

		gcd_bk = gcd(LABS(n_bra), LABS(n_ket));
		if(n_sq % gcd_bk != 0)
			*isOK = false;
		else
		{
			# solve for a and d
			startpt = n_sq % n_ket;
			i_inc=0;
			do
			{
				if((n_sq + i_inc*n_bra) % n_ket == 0)
				{
					d = i_inc;
					a = (n_sq + i_inc*n_bra)/n_ket;
					b = n_bra - a;
					c = n_ket - d;
					a2=a,b2=b,c2=c,d2=d;
					break;
				}
				else if((n_sq - i_inc*n_bra) % n_ket == 0)
				{
					d = -i_inc;
					a = (n_sq - i_inc*n_bra)/n_ket;
					b = n_bra - a;
					c = n_ket - d;
					a2=a,b2=b,c2=c,d2=d;
					break;
				}
				i_inc++;
			} 
			while(i_inc < LMAX(LABS(n_bra), LABS(n_ket)) || i_inc==1);
			
			# solve for a and d
			startpt = n_sq % n_bra;
			i_inc=0;
			do
			{
				if((i_inc*n_ket - n_sq) % n_bra == 0)
				{
					a4 = i_inc;
					d4 = (i_inc*n_ket - n_sq)/n_bra;
					b4 = n_bra - a4;
					c4 = n_ket - d4;
					break;
				}
				else if((-i_inc*n_ket - n_sq) % n_bra == 0)
				{
					a4 = -i_inc;
					d4 = (-i_inc*n_ket - n_sq)/n_bra;
					b4 = n_bra - a4;
					c4 = n_ket - d4;
					break;
				}
				i_inc++;
			} 
			while(i_inc < LMAX(LABS(n_bra), LABS(n_ket)) || i_inc==1);
			# solve for b and c
			startpt = n_sq % n_bra;
			i_inc=0;
			do
			{
				if((n_sq + i_inc*n_ket) % n_bra == 0)
				{
					b2 = i_inc;
					c2 = (n_sq + i_inc*n_bra)/n_bra;
					a2 = n_bra - b2;
					d2 = n_ket - c2;
					break;
				}
				else if((n_sq - i_inc*n_ket) % n_bra == 0)
				{
					b2 = -i_inc;
					c2 = (n_sq - i_inc*n_bra)/n_bra;
					a2 = n_bra - b2;
					d2 = n_ket - c2;
					break;
				}
				i_inc++;
			} 
			while(i_inc < LMAX(LABS(n_bra), LABS(n_ket)) || i_inc==1);

			startpt = n_sq % n_ket;
			i_inc=0;
			do
			{
				if((i_inc*n_bra - n_sq) % n_ket == 0)
				{
					c3 = i_inc;
					b3 = (i_inc*n_bra - n_sq)/n_ket;
					a3 = n_bra - b3;
					d3 = n_ket - c3;
					break;
				}
				if((n_sq + i_inc*n_bra) % n_ket == 0)
				{
					c3 = -i_inc;
					b3 = (-i_inc*n_bra - n_sq)/n_ket;
					a3 = n_bra - b3;
					d3 = n_ket - c3;
					break;
				}
				i_inc++;
			} 
			while(i_inc < LMAX(LABS(n_bra), LABS(n_ket)) || i_inc==1);


		}
		/*
		else if(true) # (n_sq >= 0)
		{
			c = myround(1.0*n_sq/n_bra);
			if((n_sq - c * n_bra) % n_ket != 0)
			{
				#*isOK = false;
				#return;
			}
			b = (c * n_bra - n_sq)/n_ket;
			d = n_ket - c;
			a = n_bra - b;

			
			a2 = myround(1.0*n_sq/n_ket);
			if((n_sq - a2 * n_ket) % n_bra != 0)
			{
				#*isOK = false;
				#return;
			}
			d2 = (a2 * n_ket - n_sq)/n_bra;
			b2 = n_bra - a2;
			c2 = n_ket - d2;
			/*
		}
		else # n_sq < 0
		{
			*/
			/*
			b3 = myround(-1.0*n_sq/n_ket);
			if((n_sq + b3 * n_ket) % n_bra != 0)
			{
				#*isOK = false;
				#return;
			}
			c3 = (b3 * n_ket + n_sq)/n_bra;
			a3 = n_bra - b3;
			d3 = n_ket - c3;

			d4 = myround(-1.0*n_sq/n_bra);
			if((n_sq + d4 * n_bra) % n_ket != 0)
			{
				#*isOK = false;
				#return;
			}
			a4 = (n_sq + d4 * n_bra) / n_ket;
			b4 = n_bra - a4;
			c4 = n_ket - d4;
		}
		*/
	}
	



	int i=0;
	vecin_arr[i++]=a;
	vecin_arr[i++]=b;
	vecin_arr[i++]=c;
	vecin_arr[i++]=d;
	vecin_arr[i++]=a2;
	vecin_arr[i++]=b2;
	vecin_arr[i++]=c2;
	vecin_arr[i++]=d2;
	vecin_arr[i++]=a3;
	vecin_arr[i++]=b3;
	vecin_arr[i++]=c3;
	vecin_arr[i++]=d3;
	vecin_arr[i++]=a4;
	vecin_arr[i++]=b4;
	vecin_arr[i++]=c4;
	vecin_arr[i++]=d4;
	/*
	vecin_arr[i++]=a;
	vecin_arr[i++]=b;
	vecin_arr[i++]=c;
	vecin_arr[i++]=d;
	vecin_arr[i++]=a2;
	vecin_arr[i++]=b2;
	vecin_arr[i++]=c2;
	vecin_arr[i++]=d2;
	*/

	return;
}

int relprime(long i, long j)
{
	if(i==1 || j == 1 || gcd(i,j)==1)
		return(true);
	else
		return(false);
}

long cgetprimeshift(long base, long pivot, long base2, long pivot2)
{

	long i, maxall = LMAX(LMAX(LABS(base), LABS(pivot)), LMAX(LABS(base2), LABS(pivot2)));

	if(gcd(LABS(base), LABS(pivot))==1  && gcd(LABS(base2), LABS(pivot2))==1)
		return 0;

	for(i=1;i<2*maxall;i++)
	{
		if(relprime(base, LABS(pivot + i)) == 1 && relprime(base, LABS(pivot - i)) == 1 && 
			relprime(base2, LABS(pivot2 - i)) == 1 && relprime(base2, LABS(pivot2 + i)) == 1)
		{
			return i;
		}
	}

	return -1;

}

void cgetfullprimeshift(long base, long pivot, long base2, long pivot2, long *i_bra, long *i_ket)
{

	long i, j, temp1, temp2, temp3, temp4, maxall = LMAX(LMAX(LABS(base), LABS(pivot)), LMAX(LABS(base2), LABS(pivot2)));

	#if(gcd(LABS(base), LABS(pivot))==1  && gcd(LABS(base2), LABS(pivot2))==1)
	#	return 0;

	for(i=0;i<maxall;i++)
	{
		temp1 = cgetprimeshift(base + i, pivot, base2 - i, pivot2);
		temp2 = cgetprimeshift(base - i, pivot, base2 + i, pivot2);
		if(temp1 > -1 && temp2 > -1)
		{
			*i_bra=temp1;
			*i_ket=0;
			return;
		}
		else
			for(j=1;j<maxall;j++)
			{
				temp1 = cgetprimeshift(base + i, pivot+j, base2 - i, pivot2-j);
				temp2 = cgetprimeshift(base - i, pivot-j, base2 + i, pivot2+j);
				temp3 = cgetprimeshift(base + i, pivot-j, base2 - i, pivot2+j);
				temp4 = cgetprimeshift(base - i, pivot-j, base2 + i, pivot2+j);
				#if(relprime(base  + i, LABS(pivot  + j)) == 1 && relprime(base  + i, LABS(pivot  - j)) == 1 && 
				#   relprime(base2 - i, LABS(pivot2 + j)) == 1 && relprime(base2 - i, LABS(pivot2 - j)) == 1)
				if(temp1>-1 && temp2>-1 && temp3>-1 && temp4>-1)
				{
					*i_bra=i;
					*i_ket=j;
					return;
				}
			}
	}

	*i_bra = *i_ket = -1;
	return;

}



void getbestampcomponents(long n_bra, long n_ket, long n_sq, double *vecin_arr, int *isOK, int freeenergy0potenergy1)
{
	int i, j, criterion;
	long a, b, c, d, abest, bbest, cbest, dbest, firsttime=true;

	getampcomponents(n_bra, n_ket, n_sq, vecin_arr, isOK);
	abest = 0;
	bbest = 0;
	cbest = 0;
	dbest = 0;

	for(i=0;i<4;i++)
	{	
		j=0;
		a = vecin_arr[4*i+j++];
		b = vecin_arr[4*i+j++];
		c = vecin_arr[4*i+j++];
		d = vecin_arr[4*i+j++];
		if((a * c - b * d == n_sq) && (a + b == n_bra) && (c + d == n_ket))
		{
			*isOK = true;
			if(firsttime==true)
			{
				abest = a;
				bbest = b;
				cbest = c;
				dbest = d;
				firsttime = false;
			}
			else
			{
				if(freeenergy0potenergy1==1)
					criterion = b*b + d*d < bbest*bbest + dbest*dbest;
				else
					criterion = a*a + b*b + c*c + d*d < abest*abest + bbest*bbest + cbest*cbest + dbest*dbest;
				if(criterion)
				{
					abest = a;
					bbest = b;
					cbest = c;
					dbest = d;
				}
			}
		}
	}
	if(*isOK==true)
	{
		int i=0;
		vecin_arr[i++]=abest;
		vecin_arr[i++]=bbest;
		vecin_arr[i++]=cbest;
		vecin_arr[i++]=dbest;
		vecin_arr[i++]=true;
	}
	else
	{
		int i=0;
		vecin_arr[i++]=0;
		vecin_arr[i++]=0;
		vecin_arr[i++]=0;
		vecin_arr[i++]=0;
		vecin_arr[i++]=false;

	}

	return;
}


void fleg(double x, double pl[], int nl)
{
	int j;
	double twox,f2,f1,d;

	pl[1]=1.0;
	pl[2]=x;
	if (nl > 2) {
		twox=2.0*x;
		f2=x;
		d=1.0;
		for (j=3;j<=nl;j++) {
			f1=d++;
			f2 += twox;
			pl[j]=(f2*pl[j-1]-f1*pl[j-2])/d;
		}
	}
}




float NormSDist(float y)
{
	double t,z,ans, x;

	x = (float) y/1.4142135623731;

	z=fabs(x);
	t=1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	ans =  (((x >= 0.0 ? -ans : -(2.0-ans))+2.)/2.);
	return (float) ans;
}






float erfcc(float x)
{
	float t,z,ans;

	z=fabs(x);
	t=1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	return (float) (x >= 0.0 ? ans : 2.0-ans);
}









Global_	  M_PI_2 1.570796326794896619231321691640 # pi/2


#http:#home.online.no/~pjacklam/notes/invnorm/ THANKS DUDE!!

Global_EPSr 1.2e-7

/*
 * Lower tail quantile for standard normal distribution function.
 *
 * This function returns an approximation of the inverse cumulative
 * standard normal distribution function.  I.e., given P, it returns
 * an approximation to the X satisfying P = Pr{Z <= X} where Z is a
 * random variable from the standard normal distribution.
 *
 * The algorithm uses a minimax approximation by rational functions
 * and the result has a relative error whose absolute value is less
 * than 1.15e-9.
 *
 * Author:      Peter J. Acklam
 * Time-stamp:  2002-06-09 18:45:44 +0200
 * E-mail:      jacklam@math.uio.no
 * WWW URL:     http:#www.math.uio.no/~jacklam
 *
 * C implementation adapted from Peter's Perl version
 */

#include <math.h>
#include <errno.h>

/* Coefficients in rational approximations. */
static const double a[] =
{
	-3.969683028665376e+01,
	 2.209460984245205e+02,
	-2.759285104469687e+02,
	 1.383577518672690e+02,
	-3.066479806614716e+01,
	 2.506628277459239e+00
};

static const double b[] =
{
	-5.447609879822406e+01,
	 1.615858368580409e+02,
	-1.556989798598866e+02,
	 6.680131188771972e+01,
	-1.328068155288572e+01
};

static const double c[] =
{
	-7.784894002430293e-03,
	-3.223964580411365e-01,
	-2.400758277161838e+00,
	-2.549732539343734e+00,
	 4.374664141464968e+00,
	 2.938163982698783e+00
};

static const double d[] =
{
	7.784695709041462e-03,
	3.224671290700398e-01,
	2.445134137142996e+00,
	3.754408661907416e+00
};

Global_LOW 0.02425
Global_HIGH 0.97575

double normsinv(double p)
{
	double q, r;

	errno = 0;

	if (p < 0+EPSr)
		return(-5e12);
	else if(p > 1-EPSr)
		return(-5e12);
	else if (p < LOW)
	{
		/* Rational approximation for lower region */
		q = sqrt(-2*log(p));
		return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else if (p > HIGH)
	{
		/* Rational approximation for upper region */
		q  = sqrt(-2*log(1-p));
		return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else
	{
		/* Rational approximation for central region */
    		q = p - 0.5;
    		r = q*q;
		return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
			(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
	}
}





#from http:#www.stat.vt.edu/~sundar/java/code/qt_js.html

Global_	  M_PI_2 1.570796326794896619231321691640 # pi/2
double qt(float p,float ndf,float lower_tail) 
{
  # Algorithm 396: Student's t-quantiles by
  # G.W. Hill CACM 13(10), 619-620, October 1970

	double eps, neg, P, q, a, b, c, d, y, x, prob;
	  if(p<=0.0 || p>=1.0 || ndf<1.0) return -1.0;
	  eps=1e-12;

	  if((lower_tail && p > 0.5) || (!lower_tail && p < 0.5)) 
	  {
		 neg = false;
		 P = 2 * (lower_tail ? (1 - p) : p);
	   }
	   else 
	   {
		 neg = true;
		 P = 2 * (lower_tail ? p : (1 - p));
	   }

	   if(abs(ndf - 2.0) < eps) 
	   {   /* df ~= 2 */
		 q=sqrt(2 / (P * (2 - P)) - 2);
	   }
	   else if (ndf < 1 + eps) 
	   {   /* df ~= 1 */
		 prob = P * M_PI_2;
		 q = cos(prob)/sin(prob);
	   }
	   else 
	   {      /*-- usual case;  including, e.g.,  df = 1.1 */
		 a = 1. / (ndf - 0.5);
		 b = 48 / (a * a);
		 c = ((20700 * a / b - 98) * a - 16) * a + 96.36;
		 d = ((94.5 / (b + c) - 3) / b + 1) * sqrt(a * M_PI_2) * ndf;
		 y = pow((double) d * P, (double) 2 / ndf);
		 if (y > 0.05 + a) {
		   /* Asymptotic inverse expansion about normal */
		   x = normsinv(0.5 * P);
		   y = x * x;
		   if (ndf < 5)
			 c += 0.3 * (ndf - 4.5) * (x + 0.6);
		   c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c;
		   y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1) * x;
		   y = a * y * y;
		   if (y > 0.002)
			 y = exp(y) - 1.0;
		   else { /* Taylor of    e^y -1 : */
			 y = (0.5 * y + 1.0) * y;
		   }
		 }
		 else {
		   y = ((1 / (((ndf + 6) / (ndf * y) - 0.089 * d - 0.822)
			   * (ndf + 2.0) * 3.0) + 0.5 / (ndf + 4.0))
			   * y - 1.0) * (ndf + 1.0) / (ndf + 2.0) + 1 / y;
		 }
		 q = sqrt(ndf * y);
	   }
	   if(neg) q = -q;
	   return q;
}




float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

#include <math.h>

float betai(float a, float b, float x)
{
	float betacf(float a, float b, float x);
	float gammln(float xx);
	#void nrerror(char error_text[]);
	float bt;

	if (x < 0.0 || x > 1.0) nrerror2("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}



#undef EPS
Global_MAXIT 100
Global_EPS 3.0e-7
Global_FPMIN 1.0e-30

float betacf(float a, float b, float x)
{
	#void nrerror(char error_text[]);
	int m,m2;
	float aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (m > MAXIT) nrerror2("a or b too big, or MAXIT too small in betacf");
	return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN



double tdist(float t, float nu)
{
	return(1.0 - betai(nu/2.0, 0.5, nu/(nu+t*t)));
}





double normtostudent(float x, float ndf)
{
	#return(x<=0.5 ? -qt(x*2.0,ndf,true) : qt((1.0-x)*2.0,ndf,true));
	return(x<=0.5 ? -qt(x,ndf,true) : qt((1.0-x),ndf,true));
}


















double gasdev(long *idum)
{
	double ran2(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if (*idum < 0) iset=0;
	if  (iset == 0) {
		do {
			v1=2.0*ran2(idum)-1.0;
			v2=2.0*ran2(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return (double) v2*fac;
	} else {
		iset=0;
		return (double) gset;
	}
}










void pre_main(int argc, char *argv[])
{

	char   str1[MAXLEN], str2[MAXLEN], str3[MAXLEN], str4[MAXLEN], str5[MAXLEN], currency[5];
	FILE   *fpin, *fpout;
	long    i, j, k, Nsteps, Nsimul, i1, i2, i3, i4, i_bk, i_dir, i_t, t_init, j1, j2, j3, j4, j_bk, j_dir, j_t;
	double temp, temp1, temp2, temp3, temp4, temp5, temp6;
	long   seed, seed_InitDeflts, sumnonzeroarcs, t_simlen;

	verbose = true;	
	std::complex <double> ctemp=0.0;






	if((fpin = fopen(argv[1], "r"))==NULL)
	{
		#if (_DEBUG) 
			printf("No input file");
		return;
	}
	/* this is how the input file looks (except for this line): */
	/*
	XYZcoord	25	25	25
	whichdirection	1		
	Nsteps	10	

	*/




	Nsteps = 100;
	seed = seed_InitDeflts = -5;
	ran2(&seed);

	t_simlen = 100000000;
	vacuumclass vacuum(NDIM, Nsteps, UNIVERSELENGTH, UNIVERSELENGTH, UNIVERSELENGTH_z);

	vacuum.grid[0][2][5][1][0]=1;
	#remember to conjugate all the stuff in the conjugate grid from the get-go then you don't have to conjugate later.
	vacuum.grid[1][2][5][1][0]=1;
	#vacuum.grid[2][0][1][1][1]=1;

	#vacuum.grid[0][1][3][0][0]=-1;
	#vacuum.grid[1][1][3][0][0]=-1;

Global_barrier_level 3
	#MAKE THE BARRIER AT LEAST 2 (or 4, or 6 or 8) POINTS THICK, TO AVOID PARITY PROBLEMS -- if it's odd, the barier is permeable from one direction.
	for(i1=0; i1<UNIVERSELENGTH; i1++)
		for(i2=barrier_level; i2< barrier_level + 0 + 0*UNIVERSELENGTH; i2++)
			for(i3=0; i3<UNIVERSELENGTH_z; i3++)
			{
				vacuum.barrier[i1][i2][i3] = 1.0;
			}



	t_init = 1;		# note we set t_init to 1 since parity would be even since we've set (i,j,k) to (25,25,25)

	if((fpout = fopen(argv[2], "w"))==NULL)
	{
		#if (_DEBUG) 
			printf("No output file");
		return;
	}
	fprintf(fpout,"\tsumampsq.real\tsumampsq.imag\tbrasumamp.real\tbrasumamp.imag\tketsumamp.real\tketsumamp.imag\tnsqsumamp.real\tnsqsumamp.imag\tbrasumabs.real\tbrasumabs.imag\tketsumabs.real\tketsumabs.imag\tnsqsumabs.real\tnsqsumabs.imag\tbra111\tket111\n"); 

	/*
	fprintf(fpout,"%d\tsumampsq%f\t%f\t\tbrasumamp%f\t%f\t\tketsumamp%f\t%f\t\t", (int) vacuum.t, 
			fprintf(fpout,"%d\tsumampsq%f\t%f\t\tbrasumamp%f\t%f\t\tketsumamp%f\t%f\t\t", (int) vacuum.t, 
				(float) vacuum.sumampsq.real(), (float) vacuum.sumampsq.imag(),
				(float) vacuum.brasumamp.real(), (float) vacuum.brasumamp.imag(),
				(float) vacuum.ketsumamp.real(), (float) vacuum.ketsumamp.imag() 
				(float) vacuum.brasumabs.real(), (float) vacuum.brasumabs.imag(),
				(float) vacuum.ketsumabs.real(), (float) vacuum.ketsumabs.imag());
	*/
	vacuum.measureamps();

	temp = vacuum.sumampsq.real();
	temp = vacuum.sumampsq.imag();
	temp = vacuum.brasumamp.real();
	temp = vacuum.brasumamp.imag();
	temp = vacuum.ketsumamp.real();
	temp = vacuum.ketsumamp.imag();
	temp = vacuum.nsqsumamp.real();
	temp = vacuum.nsqsumamp.imag();
	temp = vacuum.brasumabs.real();
	temp = vacuum.brasumabs.imag();
	temp = vacuum.ketsumabs.real();
	temp = vacuum.ketsumabs.imag();
	temp = vacuum.nsqsumabs.real();
	temp = vacuum.nsqsumabs.imag();

	#update the white squares first, then the black
	for(i_t=t_init; i_t<t_simlen; i_t++)
	{
		#once we constrain sumsq it makes sense to do both bra and ket at a given point at the same time
		#for(i_bk=0; i_bk<2; i_bk++)
			for(i1=0; i1<UNIVERSELENGTH; i1++)
				for(i2=0; i2<UNIVERSELENGTH; i2++)
					for(i3=0; i3<UNIVERSELENGTH_z; i3++)
					{
						#if((i1>=25 && i1<=26) && (i2>=25 && i2<=25)  && (i3>=25 && i3<=25))
						#	i1 = i1;
						if((i1+i2+i3+i_t)%2==0) # parity positive
							continue;
						else
						{
							/*
							if(i1==2 && i2==1 && i3==UNIVERSELENGTH-1)
								temp = temp;
							if(i1==1 && i2==0 && i3==0 && i_t==6)
								temp = temp;
							if(i1==2 && i2==1 && i3==0 && (i_t==4 || i_t==5))
								temp = temp;
							*/
							if(vacuum.isempty_pt(0, i1, i2, i3, BACKSHIFT) == false ) 
							{
								temp1 = vacuum.realsum_pt(0,i1,i2,i3,BACKSHIFT);
								#temp2 = vacuum.realsum_pt(1,i1,i2,i3,ZEROSHIFT);
								#temp1 = vacuum.realsum_pt(0,i1,i2,i3,BACKSHIFT);
								#temp2 = vacuum.realsum_pt(1,i1,i2,i3,BACKSHIFT);
								
								#temp3 = temp1 * temp2;
								vacuum.generateoutflow_from(0, i1, i2, i3, &seed);

								temp4 = vacuum.realsum_pt(0,i1,i2,i3,ZEROSHIFT);
								#temp5 = vacuum.realsum_pt(1,i1,i2,i3,ZEROSHIFT);
								#temp4 = vacuum.realsum_pt(0,i1,i2,i3,BACKSHIFT);
								#temp5 = vacuum.realsum_pt(1,i1,i2,i3,BACKSHIFT);
								

								if(fabs(temp1 - temp4) > .000000001)
									temp1 = temp1;
							}
							if(vacuum.isempty_pt(1, i1, i2, i3, BACKSHIFT) == false)
								vacuum.generateoutflow_from(1, i1, i2, i3, &seed);
						}
					}
			#now wipe the off-parity grid points clean
			for(i1=0; i1<UNIVERSELENGTH; i1++)
				for(i2=0; i2<UNIVERSELENGTH; i2++)
					for(i3=0; i3<UNIVERSELENGTH_z; i3++)
					{
						#if((i1>=25 && i1<=26) && (i2>=25 && i2<=25)  && (i3>=25 && i3<=25))
						#	i1 = i1;
						if((i1+i2+i3+i_t)%2==1) # parity negative
							continue;
						else
						{
							/*
							if(i1==1 && i2==1 && i3==1 && i_t==3)
								temp = temp;
							if(i1==1 && i2==1 && i3==0*UNIVERSELENGTH - -1)
								i1=i1;
							*/
							for(i_dir=0;i_dir<NDIM2;i_dir++)
							{
								/*
								ctemp = vacuum.grid[0][i_dir][i1][i2][i3];
								temp1 = ctemp.real();
								ctemp = vacuum.grid[1][i_dir][i1][i2][i3].real();
								temp2 = ctemp.real();
								temp3 = vacuum.sumsq_pt(i1, i2, i3, ZEROSHIFT);
								*/

								vacuum.grid[0][i_dir][i1][i2][i3] = 0.0;
								vacuum.grid[1][i_dir][i1][i2][i3] = 0.0;

								/*
								ctemp = vacuum.outflow(0, i_dir, i1, i2, i3);
								temp4 = ctemp.real();
								ctemp = vacuum.outflow(1, i_dir, i1, i2, i3);
								temp5 = ctemp.real();
								temp6 = vacuum.outflow_sumsq(i_dir, i1, i2, i3);
								if(fabs(temp1 - temp4) > .000000001)
									temp1 = temp1;
								*/
							}
						}
					}
					
					
		vacuum.t++;
		vacuum.measureamps();
		temp = vacuum.sumampsq.real();
		temp = vacuum.sumampsq.imag();
		temp = vacuum.brasumamp.real();
		temp = vacuum.brasumamp.imag();
		temp = vacuum.ketsumamp.real();
		temp = vacuum.ketsumamp.imag();
		temp = vacuum.nsqsumamp.real();
		temp = vacuum.nsqsumamp.imag();

		temp = vacuum.brasumabs.real();
		temp = vacuum.brasumabs.imag();
		temp = vacuum.ketsumabs.real();
		temp = vacuum.ketsumabs.imag();
		temp = vacuum.nsqsumabs.real();
		temp = vacuum.nsqsumabs.imag();
	
		if(i_t%1==0)
		{
			fprintf(fpout,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t\t\n", (int) vacuum.t, 
				(float) vacuum.sumampsq.real(), (float) vacuum.sumampsq.imag(),
				(float) vacuum.brasumamp.real(), (float) vacuum.brasumamp.imag(),
				(float) vacuum.ketsumamp.real(), (float) vacuum.ketsumamp.imag(), 
				(float) vacuum.nsqsumamp.real(), (float) vacuum.nsqsumamp.imag(), 
				(float) vacuum.brasumabs.real(), (float) vacuum.brasumabs.imag(),
				(float) vacuum.ketsumabs.real(), (float) vacuum.ketsumabs.imag(),
				(float) vacuum.nsqsumabs.real(), (float) vacuum.nsqsumabs.imag(), (float) vacuum.grid[0][0][1][1][1].real(), vacuum.grid[1][0][1][1][1].real() );
			if(((long) vacuum.t) % 10 == 0)
				printf("%d\n", (long) vacuum.t);
		}
		
		
		for(i1=0;i1<UNIVERSELENGTH_z;i1++)
		{
			fprintf(fpout, "\nz=%d\tt=%d\n", i1, i_t);
			vacuum.printamps(fpout,3,i1,0,0);
		}		
		

		fflush(fpout);
	}


}




int _tmain(int argc, _TCHAR* argv[])
{

	pre_main(argc, argv);
	long idum=-5;


	gasdev(&idum);
	return 0;
}




/*======================================================================

Author: Jin Chu Wu  Ph.D.
        National Institute of Standards and Technology

The class map,

  GetData

  Util

                                  SRS   StatPackage
                                   ^     ^
                                   |     |
              TwoArrays          BootstrapCI
                  ^                   ^
                  |                   |
                  ---------------------
                           |
                           |
                       CellOperations
                           ^
                           |
                          CIS

The options of the command line are,

         [-b]    : Cell Image Segmentation -- bootstrap method
         [-a]    : Cell Image Segmentation -- analytical method
         [-c]    : Compare MERs of two CIS algorithms
         [-h]    : Help flag

======================================================================*/

#include <stdio.h>
// stdlib.h contains definitions of EXIT_SUCCESS and EXIT_FAILURE.
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <malloc.h>
#include <math.h>
#include <unistd.h>
#include <dirent.h>
#include <ctype.h>
#include <limits.h>

typedef char Boolean;
#ifndef TRUE
#define TRUE (Boolean)1
#endif
#ifndef FALSE
#define FALSE (Boolean)0
#endif

#define STRINGSIZE                256
#define bsSigLevel               0.05
#define ZSCORE               1.959964
#define R_PATH "R"

#define max2(a,b)	((a) > (b) ? (a) : (b))
#define min2(a,b)	((a) > (b) ? (b) : (a))

int intCompDescend (const void * a, const void * b)
{
  int * x, * y;

  x = (int *) a;
  y = (int *) b;

  if ((* x) > (* y))
    return -1;
  else if ((* x) == (* y))
    return 0;
  else
    return 1;
}

int intCompAscend (const void * a, const void * b)
{
  int * x, * y;

  x = (int *) a;
  y = (int *) b;

  if ((* x) > (* y))
    return 1;
  else if ((* x) == (* y))
    return 0;
  else
    return -1;
}

int doubleCompAscend (const void * a, const void * b)
{
  double * x, * y;

  x = (double *) a;
  y = (double *) b;

  if ((* x) > (* y))
    return 1;
  else if ((* x) == (* y))
    return 0;
  else
    return -1;
}

// ========== class GetData ==========

class GetData {

  char CISFlag, *inputCISFilename, *inputCISFilename2;
  int whichRate, bootstrapReplicationNum, inputThresholdCIS, iterationNum;
  double confidenceLevel;
  long rngSeed;

  int nG, nA, ng, na;
  double st, se, mu, st1, se1, st2, se2, corCoef;

 public:

  char flag;

  GetData (int, char **);

  friend class Util;
  friend class SRS;
  friend class TwoArrays;
  friend class BootstrapCI;
  friend class CellOperations;
  friend class CIS;

};

GetData::GetData (int argc, char ** argv)
{
  rngSeed = 0;
  st = se = mu = st1 = se1 = st2 = se2 = corCoef = 0.0;
  nG = nA = ng = na = 0;

  if (strncmp (argv [1], "-b", 2) == 0)
    {
      flag = 'b';

      if ((argc == 5) || (argc == 6))
	{
	  // Multiple cells

	  CISFlag = 'n';

	  iterationNum = atoi (argv [3]);

	  inputCISFilename = argv [4];

	  if (argc == 5)
	    rngSeed = 0;
	  else
	    rngSeed = atol (argv [5]);
	}
      else if ((argc == 7) || (argc == 8))
	{
	  // One cell

	  CISFlag = '1';

	  iterationNum = 0;

	  nG = atoi (argv [3]);
	  nA = atoi (argv [4]);

	  ng = atoi (argv [5]);
	  na = atoi (argv [6]);

	  if (argc == 7)
	    rngSeed = 0;
	  else
	    rngSeed = atol (argv [7]);
	}
      else
	{
	  printf ("The command line: cellImageSeg -b whichRate nG nA ng na OR \
cellImageSeg -b whichRate nG nA ng na rngSeed OR \
cellImageSeg -b whichRate inputFileName OR \
cellImageSeg -b whichRate inputFileName rngSeed\n");

	  exit (1);
	}
    }
  else if (strncmp (argv [1], "-a", 2) == 0)
    {
      flag = 'a';

      if (argc == 4)
	{
	  inputCISFilename = argv [3];
	}
      else
	{
	  printf ("The command line: \
cellImageSeg -a whichRate inputFileName\n");

	  exit (1);
	}
    }
  else if (strncmp (argv [1], "-c", 2) == 0)
    {
      flag = 'c';

      if (argc == 5)
	{
	  inputCISFilename = argv [3];
	  inputCISFilename2 = argv [4];
	}
      else
	{
	  printf ("The command line: \
cellImageSeg -c whichRate inputFileName1 inputFileName2\n");

	  exit (1);
	}
    }
  else if (strncmp (argv [1], "-h", 2) == 0)
    {
printf ("         [-b]    : Cell Image Segmentation -- bootstrap method\n");
printf ("         [-a]    : Cell Image Segmentation -- analytical method\n");
      printf ("         [-h]    : Help flag\n");
      exit (1);
    }
  else
    {
      printf ("There is something wrong with the flag.\n");
      exit (1);
    }

  whichRate = atoi (argv [2]);

  if ((flag == 'a') && (whichRate != 1))
    {
    printf ("If using the analytical method, the average MER must be used.\n");

      exit (1);
    }

  bootstrapReplicationNum = 2000;
  confidenceLevel = 0.95;
  inputThresholdCIS = 1;
}

// ========== class Util ==========

class Util {

  const GetData * data;

  static void callSystem (const char *, ...);
  template<class R> static void allocateOneArray (int, R **);
  template<class R> static void bzeroOneArray (int, R *);
  template<class R, class S> static void allocateTwoArrays (int, R **,
							    int, S **);
  template<class R, class S> static void bzeroTwoArrays (int, R *, int, S *);

 public:

  Util (const GetData& inputData)
  {
    data = &inputData;
  }

  friend class TwoArrays;
  friend class BootstrapCI;
  friend class CellOperations;
  friend class CIS;

};

void Util::callSystem (const char * counter, ...)
{
  char * term, systemString [STRINGSIZE];
  va_list ap;

  bzero ((char *) systemString, STRINGSIZE * sizeof (char));

  va_start (ap, counter);

  while (*counter++)
    {
      term = va_arg (ap, char *);

      strcat (systemString, term);
    }

  va_end (ap);

  system (systemString);
}

template<class R> void Util::allocateOneArray (int num, R **array)
{
  R *tmpArray;

  tmpArray = (R *) malloc (num * sizeof (R));

  *array =  tmpArray;
}

template<class R> void Util::bzeroOneArray (int num, R *array)
{
  bzero ((R *) array, num * sizeof (R));
}

template<class R, class S> void Util::allocateTwoArrays
(int num1, R ** array1, int num2, S ** array2)
{
  R * tmpArray1;
  S * tmpArray2;

  tmpArray1 = (R *) malloc (num1 * sizeof (R));
  tmpArray2 = (S *) malloc (num2 * sizeof (S));

  * array1 =  tmpArray1;
  * array2 =  tmpArray2;
}

template<class R, class S> void Util::bzeroTwoArrays
(int num1, R * array1, int num2, S * array2)
{
  bzero ((R *) array1, num1 * sizeof (R));
  bzero ((S *) array2, num2 * sizeof (S));
}

// ========== class SRS ==========

class SRS {

  const GetData * data;
  long RNGSEED;

  double ran2(long *idum);
  double ranDoubleSRS ();
  Boolean ranBooleanSRS (double);
  int ranIntSRS (int);
  long createRNGSEED ();

 protected:

  void obtainRNGSEED ();
  template<class S, class T> void arrayToArrayWRSRS (int, int, S *, T *);
  template<class R, class S> void synchronizedWRSRS (int, R *, R *, S *, S *);
  void numberToArrayWORSRS (int, int, int *, char *);
  template<class S, class T> void arrayToArrayWORSRS (int, int, S *, T *,
						      char *);

 public:

  SRS (const GetData& inputData)
  {
    data = &inputData;
  }

};

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

// Change "ran2" and "temp" from "float" to "double".

double SRS::ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

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
	else return temp;
}

double SRS::ranDoubleSRS ()
{
  return (ran2 (&RNGSEED));
}

Boolean SRS::ranBooleanSRS (double prob)
{
  return ((Boolean) (ranDoubleSRS () < prob));
}

int SRS::ranIntSRS (int LENGTHSRS)
{
  return ((int) (ranDoubleSRS () * (double) LENGTHSRS));
}

void SRS::obtainRNGSEED ()
{
  if (data->rngSeed > 0)
    {
      RNGSEED = - (data->rngSeed);
    }
  else
    {
      RNGSEED = createRNGSEED ();
    }

  printf ("RNGSEED = %ld\n", - RNGSEED);
}

long SRS::createRNGSEED ()
{
  unsigned long int bigNumber = 0;
  FILE * fp;

  while (1)
    {
      if ((fp = popen ("ps -elf | cksum | awk '{print $1}'", "r")) == NULL)
	{
	  printf ("It cannot read output from ps command");
	  exit (1);
	}
      fscanf (fp, "%lu", &bigNumber);
      pclose (fp);

      // The following condition is made based on the function ran2.

      if (bigNumber <= 4294967295)
	{
	  if ((__WORDSIZE != 64) && (bigNumber > 2147483647))
	    bigNumber /= 2;

	  return (- bigNumber);
	}
    }
}

template<class S, class T> void SRS::arrayToArrayWRSRS
(int dataNum, int sampleNum, S * data, T * samples)
{
  register int i, index;

  for (i = 0; i < sampleNum; i++)
    {
      while (1)
	{
	  index = ranIntSRS (dataNum);

	  if (index != dataNum)
	    {
	      samples [i] = data [index];

	      break;
	    }
	}
    }
}

template<class R, class S> void SRS::synchronizedWRSRS
(int num, R * rawData1, R * sample1, S * rawData2, S * sample2)
{
  register int count, index;

  count = 0;

  while (1)
    {
      index = ranIntSRS (num);

      if (index != num)
	{
	  sample1 [count] = rawData1 [index];
	  sample2 [count] = rawData2 [index];

	  if (++count == num)
	    break;
	}
    }
}

void SRS::numberToArrayWORSRS (int numData, int numSample,
			       int *sampleArray, char *flag)
{
  register int count, index;

  count = 0;

  while (1)
    {
      index = ranIntSRS (numData);

      if (index != numData)
	{
	  if (flag [index] == FALSE)
	    {
	      flag [index] = TRUE;

	      sampleArray [count] = index;

	      if (++count == numSample)
		break;
	    }
	}
    }
}

template<class S, class T> void SRS::arrayToArrayWORSRS
(int numData, int numSample, S *data, T *samples, char *srsFlag)
{
  register int count, index;

  count = 0;

  while (1)
    {
      index = ranIntSRS (numData);

      if (index != numData)
	{
	  if (srsFlag [index] == FALSE)
	    {
	      srsFlag [index] = TRUE;

	      samples [count] = data [index];

	      if (++count == numSample)
		break;
	    }
	}
    }
}

// ========== class StatPackage ==========

class StatPackage {

  static double computeNumeratorStDev (int, double *);
  double runRNorm (double, const char *);

 protected:

  void runRPercentile (const char *, double, double,
		       double *, double *, double *);
  static double computeMean (int, double *);
  double computeUnbiasedStDev (int, double *);
  double computeCorrelationCoefficient (int, double *, double *);
  double computeCV (int, double *);
  double probToZScore (double);
  double ZScoreToProb (double);
  double computePercentile (int, double *, double, char);
  void OneSample_ZTest (double, double, double, double *, double *);
  void TwoSamples_ZTest (double, double, double, double, double,
			 double *, double *);

 public:

  StatPackage ()
  {
  }

};

double StatPackage::computeMean (int num, double * dataSet)
{
  int i;
  double sum;

  sum = 0.0;

  for (i = 0; i < num; i++)
    sum += dataSet [i];

  return (sum / (double) num);
}

double StatPackage::computeNumeratorStDev (int num, double * dataSet)
{
  int i;
  double mean, numerator, dif;

  mean = computeMean (num, dataSet);

  numerator = 0.0;

  for (i = 0; i < num; i++)
    {
      dif = dataSet [i] - mean;
      numerator += dif * dif;
    }

  return (numerator);
}

double StatPackage::computeUnbiasedStDev (int num, double * dataSet)
{
  return (sqrt (computeNumeratorStDev (num, dataSet) / (double) (num - 1)));
}

void StatPackage::runRPercentile (const char *fileName,
				  double alphaLo, double alphaUp,
				double *ciLo2, double *medium, double *ciUp2)
{
  const char *RCommand1, *RCommand2, *ct0 = " \n\t[]";
  char buffer [STRINGSIZE], *tmp;
  FILE *fpRunR;

  // Write an R script.

  RCommand1 = "quantile (a[,1], probs = c(";
  RCommand2 = "), type = ";

  fpRunR = fopen ("__RUN_TEMP_R_SCRIPT__", "w");
  fprintf (fpRunR, "%s\n\n", "#!/bin/csh -f");
  fprintf (fpRunR, "%s%s\n", R_PATH, " --slave --no-save <<\\\\");
  fprintf (fpRunR, "%s", "a <- read.table(");
  fprintf (fpRunR, "\"%s\"", fileName);
  fprintf (fpRunR, "%s\n", ")");

  fprintf (fpRunR, "%s", RCommand1);
  fprintf (fpRunR, "%5.3f, %5.3f, %5.3f", alphaLo, 0.5, alphaUp);
  fprintf (fpRunR, "%s", RCommand2);
  fprintf (fpRunR, "%i)\n", 2);

  fprintf (fpRunR, "%s\n%s\n", "q()", "\\\\");
  fclose (fpRunR);

  // Run an R script.

  system ("chmod 755 __RUN_TEMP_R_SCRIPT__");
  system ("__RUN_TEMP_R_SCRIPT__ > __TEMP_R_RESULTS__");
  system ("rm -f __RUN_TEMP_R_SCRIPT__");

  // Read the R output.
  // The format is
  //      2.5%       50%     97.5%
  // 0.9928727 0.5000000 0.9938568

  fpRunR = fopen ("__TEMP_R_RESULTS__", "r");
  bzero ((char *) buffer, STRINGSIZE * sizeof (char));

  fgets (buffer, STRINGSIZE, fpRunR);
  fgets (buffer, STRINGSIZE, fpRunR);
  tmp = NULL;
  tmp = (char *) strtok (buffer, ct0);
  *ciLo2 = atof (tmp);
  tmp = (char *) strtok (NULL, ct0);
  *medium = atof (tmp);
  tmp = (char *) strtok (NULL, ct0);
  *ciUp2 = atof (tmp);

  fclose (fpRunR);
  system ("rm -f __TEMP_R_RESULTS__");
}

double StatPackage::computeCorrelationCoefficient (int dataNum,
						   double * dataSet1,
						   double * dataSet2)
{
  int i;
  double mean1, mean2, sum, numeratorStDev1, numeratorStDev2, r;

  mean1 = computeMean (dataNum, dataSet1);
  mean2 = computeMean (dataNum, dataSet2);

  sum = 0.0;

  for (i = 0; i < dataNum; i++)
    sum += (dataSet1 [i] - mean1) * (dataSet2 [i] - mean2);

  numeratorStDev1 = computeNumeratorStDev (dataNum, dataSet1);
  numeratorStDev2 = computeNumeratorStDev (dataNum, dataSet2);

  r = sum / (sqrt (numeratorStDev1) * sqrt (numeratorStDev2));

  return (r);
}

double StatPackage::computeCV (int num, double * seq)
{
  return (computeUnbiasedStDev (num, seq) / computeMean (num, seq));
}

double StatPackage::probToZScore (double prob)
{
  return (runRNorm (prob, "qnorm (a, 0, 1)"));
}

double StatPackage::ZScoreToProb (double zscore)
{
  return (runRNorm (zscore, "pnorm (a, 0, 1)"));
}

double StatPackage::runRNorm (double d, const char * RCommand)
{
  const char *ct0 = " \n\t[]";
  char buffer [STRINGSIZE], *tmp;
  FILE *fpRunR;

  // Write an R script.

  fpRunR = fopen ("__RUN_TEMP_R_SCRIPT__", "w");
  fprintf (fpRunR, "%s\n\n", "#!/bin/csh -f");
  fprintf (fpRunR, "%s%s\n", R_PATH, " --slave --no-save <<\\\\");
  fprintf (fpRunR, "%s", "a <- c(");
  fprintf (fpRunR, " %f ", d);
  fprintf (fpRunR, "%s\n", ")");
  fprintf (fpRunR, "%s\n", RCommand);
  fprintf (fpRunR, "%s\n%s\n", "q()", "\\\\");
  fclose (fpRunR);

  // Run an R script.

  system ("chmod 755 __RUN_TEMP_R_SCRIPT__");
  system ("__RUN_TEMP_R_SCRIPT__ > __TEMP_R_RESULTS__");
  system ("rm -f __RUN_TEMP_R_SCRIPT__");

  // Read the R output.
  // The format is
  // [1] value

  fpRunR = fopen ("__TEMP_R_RESULTS__", "r");
  bzero ((char *) buffer, STRINGSIZE * sizeof (char));
  fgets (buffer, STRINGSIZE, fpRunR);
  tmp = NULL;
  tmp = (char *) strtok (buffer, ct0);
  tmp = (char *) strtok (NULL, ct0);
  fclose (fpRunR);
  system ("rm -f __TEMP_R_RESULTS__");

  return (atof (tmp));
}

double StatPackage::computePercentile (int num, double * dataSet, double alpha,
				       char flag)
{
  int index;

  if ((alpha < 0.5) || ((alpha == 0.5) && (flag == '1')))
    {
      index = (int) floor (num * alpha);
      return (dataSet [index]);
/*
      index = (int) floor ((num + 1) * alpha);
      return (dataSet [index]);

      return (dataSet [index - 1]);
*/
    }
  else
    {
      index = (int) floor (num * (1.0 - alpha));
      return (dataSet [num - index - 1]);
/*
      index = (int) floor ((num + 1) * (1.0 - alpha));
      return (dataSet [num + 1 - index]);

      return (dataSet [num - index]);
*/
    }
}

void StatPackage::OneSample_ZTest (double st, double se, double mu,
				   double * ZScore, double * pValue)
{
  *ZScore = fabs (st - mu) / se;

  *pValue = 2.0 * (1.0 - ZScoreToProb (*ZScore));
}

void StatPackage::TwoSamples_ZTest (double st1, double se1,
				    double st2, double se2, double corCoef,
				    double * ZScore, double * pValue)
{
  *ZScore = fabs (st1 - st2) /
    sqrt (se1 * se1 + se2 * se2 - 2.0 * corCoef * se1 * se2);

  *pValue = 2.0 * (1.0 - ZScoreToProb (*ZScore));
}

// ========== class BootstrapCI ==========

class BootstrapCI : public StatPackage, public SRS {

  const GetData * data;
  int * bsSampleGT, * bsSampleAlg, inputThresholdCIS;
  double * MERReplications;

 protected:

  int maxArrayLength, bsReplicationNum;
  double CILo, CIUp, BCaCILo, BCaCIUp;

  void initializeArraysBS (int, int);
  double bootstrapSEOfMERate (char, int, int, int, int, int, int *, int *,
			   double dummy (int, double, double, double, double));
  void computeBCaCI (int, int, int, int, int, int *, int *, double,
		     double dummy (int, double, double, double, double));

 public:

  BootstrapCI (const GetData& inputData) : StatPackage (), SRS (inputData)
  {
    CILo = CIUp = BCaCILo = BCaCIUp = 0.0;

    data = &inputData;

    bsReplicationNum = data->bootstrapReplicationNum;
    inputThresholdCIS = data->inputThresholdCIS;
  }

};

void BootstrapCI::initializeArraysBS (int num, int bsReplicationNum)
{
  Util::allocateOneArray (num, &bsSampleGT);
  Util::allocateOneArray (num, &bsSampleAlg);
  Util::allocateOneArray (bsReplicationNum, &MERReplications);
}

double BootstrapCI::bootstrapSEOfMERate (char flag, int whichRate,
					 int NG, int NA, int Ng, int Na,
					 int *GT, int *ALG,
			    double dummy (int, double, double, double, double))
{
  int i, j, count;
  double dg, da, dummyMedium;
  FILE *fp, *fpR;

  Util::bzeroOneArray (maxArrayLength, bsSampleGT);
  Util::bzeroOneArray (maxArrayLength, bsSampleAlg);
  Util::bzeroOneArray (bsReplicationNum, MERReplications);

  if (flag == '1')
    {
      fp = fopen ("bootstrapReps_CISMER", "w");
      fpR = fopen ("__REPLICATIONS_TEMP_FILE_1__", "w");
    }

  for (i = 0; i < bsReplicationNum; i++)
    {
      // One-sample bootstrap on algorithm-detected cell only

      if ((Ng == NG) && (Na == NA))
	{
	  return (0.0);
	}
      else if ((Ng == 0) && (Na == 0))
	{
	  return (0.0);
	}
      else if ((Ng != 0) && (Na == 0))
	{
	  while (1)
	    {
	      arrayToArrayWRSRS (NG, NG, GT, bsSampleGT);

	      count = 0;
	      for (j = 0; j < NG; j++)
		{
		  if (bsSampleGT [j] > inputThresholdCIS)
		    count++;
		}

	      if (count <= NA)
		{
		  dg = (double) (NG - count);
		  da = (double) (NA - count);

		  break;
		}
	    } // while (1)
	}
      else
	{
	  while (1)
	    {
	      arrayToArrayWRSRS (NA, NA, ALG, bsSampleAlg);

	      count = 0;
	      for (j = 0; j < NA; j++)
		{
		  if (bsSampleAlg [j] < inputThresholdCIS)
		    count++;
		}

	      if (count <= NG)
		{
		  dg = (double) (NG - count);
		  da = (double) (NA - count);

		  break;
		}
	    } // while (1)
	  } // if

      MERReplications [i] = dummy (whichRate,
				   (double) NG, (double) NA, dg, da);

      if (flag == '1')
	{
	  fprintf (fp, "%0.18f\n", MERReplications [i]);
	  fprintf (fpR, "%0.18f\n", MERReplications [i]);
	}
    } // for (i = 0; i < bsReplicationNum; i++)

  if (flag == '1')
    {
      fclose (fp);
      fclose (fpR);

      runRPercentile ("__REPLICATIONS_TEMP_FILE_1__",
		      bsSigLevel / 2.0, 1.0 - bsSigLevel / 2.0,
		      &CILo, &dummyMedium, &CIUp);

      system ("rm -f __REPLICATIONS_TEMP_FILE_1__");
    }  

  return (computeUnbiasedStDev (bsReplicationNum, MERReplications));
}

void BootstrapCI::computeBCaCI (int whichRate,
				int NG, int NA, int Ng, int Na,
				int *GT, int *ALG, double plugInMER,
			    double dummy (int, double, double, double, double))
{
  int i, num;
  double z;

  // Compute the bias-correction z.

  num = 0;
  for (i = 0; i < bsReplicationNum; i++)
    {
      if (MERReplications [i] < plugInMER)
	{
	  num++;
	}
    }

  // Modify it. Otherwise, if num = 0, then z = probToZScore (0.0) = -inf.
  // This can be tested using "% cellImageSeg -a 0 258 295 1 38 1819554413".

  if (num == 0)
    {
      for (i = 0; i < bsReplicationNum; i++)
	{
	  if (MERReplications [i] <= plugInMER)
	    {
	      num++;
	    }
	}
    }

  z = probToZScore ((double) num / (double) bsReplicationNum);

  // Compute the acceleration a.

  double a, *jkMER, temp1, temp2, jkMERAver, dif, difSqu, squareMER, cubeMER;

  Util::allocateOneArray (NA, &jkMER);

  temp1 = dummy (whichRate,
		 (double) NG, (double) (NA - 1),
		 (double) Ng, (double) (Na - 1));

  temp2 = dummy (whichRate,
		 (double) NG, (double) (NA - 1),
		 (double) (Ng + 1), (double) Na);

  for (i = 0; i < NA; i++)
    {
      if (i < Na)
	{
	  jkMER [i] = temp1;
	}
      else
	{
	  jkMER [i] = temp2;
	}
    }

  jkMERAver = computeMean (NA, jkMER);
  //printf ("%0.18f\n", (Na * temp1 + (NA - Na) * temp2) / (double) NA);

  squareMER = 0.0;
  cubeMER = 0.0;

  for (i = 0; i < NA; i++)
    {
      dif = jkMERAver - jkMER [i];
      difSqu = dif * dif;

      squareMER += difSqu;
      cubeMER += difSqu * dif;
    }

  a = cubeMER / 6.0 / sqrt (squareMER * squareMER * squareMER);

  // Compute the BCa confidence interval (BCaCILo, BCaCIUp).

  double halfSigLevel, normZLo, normZUp, tmp, alphaLo, alphaUp, dummyMedium;
  FILE *fpTemp;

  halfSigLevel = bsSigLevel / 2.0;

  normZLo = probToZScore (halfSigLevel);
  normZUp = probToZScore (1.0 - halfSigLevel);

  tmp = z + normZLo;
  alphaLo = ZScoreToProb (z + tmp / (1.0 - a * tmp));
  tmp = z + normZUp;
  alphaUp = ZScoreToProb (z + tmp / (1.0 - a * tmp));

  fpTemp = fopen ("__REPLICATIONS_TEMP_FILE_2__", "w");

  for (i = 0; i < bsReplicationNum; i++)
    {
      fprintf (fpTemp, "%0.18f\n", MERReplications [i]);
    }

  runRPercentile ("__REPLICATIONS_TEMP_FILE_2__",
		  alphaLo, alphaUp, &BCaCILo, &dummyMedium, &BCaCIUp);

  system ("rm -f __REPLICATIONS_TEMP_FILE_2__");
  fclose (fpTemp);

  /*
  qsort (MERReplications, bsReplicationNum, sizeof (double), doubleCompAscend);

  BCaCILo = computePercentile (bsReplicationNum, MERReplications, alphaLo, '1');
  BCaCIUp = computePercentile (bsReplicationNum, MERReplications, alphaUp, '2');
  */

  if (BCaCILo > BCaCIUp)
    {
 printf ("The confidence interval endpoints using BCa are not appropriate.\n");
 //exit (1);
    }
}

// ========== class TwoArrays ==========

class TwoArrays {

  const GetData * data;

  char tarFarFreqFlag;
  int inputThresholdCIS;

  void assignValuesToGtAndAlg (int, int, int *, int, int, int *);

 protected:

  int tarNum, farNum, * tar, * far, * tarFreq, * farFreq;

  void createGtAndAlgArrays (int, int, int, int, int **, int **);
  void assignValuesToGT (int, int, int *);
  void assignValuesToAlg (int, int, int *);
  void compareTwoArrays (int, double *, double *);

 public:

  TwoArrays (const GetData& inputData)
  {
    data = &inputData;

    tarFarFreqFlag = FALSE;

    tarNum = farNum = 0;
    tar = far = NULL;

    inputThresholdCIS = data->inputThresholdCIS;
  }

  ~TwoArrays ()
  {
    free (tar);
    free (far);

    if (tarFarFreqFlag)
      {
	free (tarFreq);
	free (farFreq);
      }
  }

};

void TwoArrays::createGtAndAlgArrays (int nG, int nA, int ng, int na,
				      int **GT, int **Alg)
{
  Util::allocateTwoArrays (nG, GT, nA, Alg);

  assignValuesToGtAndAlg (nG, ng, *GT, nA, na, *Alg);
}

void TwoArrays::assignValuesToGtAndAlg (int nG, int ng, int *GT,
					int nA, int na, int *Alg)
{
  int i;

  Util::bzeroTwoArrays (nG, GT, nA, Alg);

  for (i = 0; i < nG; i++)
    {
      if (i < nG - ng)
	GT [i] = inputThresholdCIS + 1;
      else
	GT [i] = inputThresholdCIS - 1;
    }

  for (i = 0; i < nA; i++)
    {
      if (i < na)
	Alg [i] = inputThresholdCIS + 1;
      else
	Alg [i] = inputThresholdCIS - 1;
    }
}

void TwoArrays::assignValuesToGT (int nG, int ng, int *gt)
{
  int i;

  Util::bzeroOneArray (nG, gt);

  for (i = 0; i < nG; i++)
    {
      if (i < ng)
	gt [i] = inputThresholdCIS - 1;
      else
	gt [i] = inputThresholdCIS + 1;
    }
}

void TwoArrays::assignValuesToAlg (int nA, int na, int *Alg)
{
  int i;

  Util::bzeroOneArray (nA, Alg);

  for (i = 0; i < nA; i++)
    {
      if (i < na)
	Alg [i] = inputThresholdCIS + 1;
      else
	Alg [i] = inputThresholdCIS - 1;
    }
}

void TwoArrays::compareTwoArrays (int numCells,
				  double *array1, double *array2)
{
  int i, numS, numG, numE;
  FILE *fp;

  if ((fp = fopen ("comparisonOfMERsOfTwoAlgs", "w")) == NULL)
    {
      printf ("The output file of comparison cannot be opened.\n");
      exit (1);
    }

  numS = numG = numE = 0;

  for (i = 0; i < numCells; i++)
    {
      fprintf (fp, "%3d, %f %f", i + 1, array1 [i], array2 [i]);

      if (array1 [i] < array2 [i])
	{
	  fprintf (fp, "   <\n");
	  numS++;
	}
      else if (array1 [i] > array2 [i])
	{
	  fprintf (fp, "       >\n");
	  numG++;
	}
      else if (array1 [i] == array2 [i])
	{
	  fprintf (fp, "           =\n");
	  numE++;
	}
    }

  fprintf (fp, "The number of < = %i\n", numS);
  fprintf (fp, "The number of > = %i\n", numG);
  if (numE != 0)
    fprintf (fp, "The number of = = %i\n", numE);

  fclose (fp);
}

// ========== class CellOperations ==========

class CellOperations: public TwoArrays, public BootstrapCI {

  const GetData * data;

  int whichRate;

  static double computeWeightedRate (double, double);
  static double computeAverageRate (double, double);

 protected:

  typedef struct {
    int nG, ng, nA, na, *gt, *alg;
  } CELL;

  #define nG(c,i)       (((c)+(i))->nG)
  #define ng(c,i)       (((c)+(i))->ng)
  #define nA(c,i)       (((c)+(i))->nA)
  #define na(c,i)       (((c)+(i))->na)
  #define gt(c,i)       (((c)+(i))->gt)
  #define gtP(c,i,j)    (((c)+(i))->gt+(j))
  #define gtI(c,i,j)    (((c)+(i))->gt[j])
  #define alg(c,i)      (((c)+(i))->alg)
  #define algP(c,i,j)   (((c)+(i))->alg+(j))
  #define algI(c,i,j)   (((c)+(i))->alg[j])

  CELL *cells;
  char *inputCISFilename, *inputCISFilename2;
  int numCells;
  double totalProb, tpSE;

  static double computeMERate (int, double, double, double, double);

  int countNumCells (char *);
  void buildCellArrays (int, char *, CELL *);
  int determineMaxArrayLength (int, CELL *);
  void computeTotalProbAndSE (char, int, CELL *, double *, double *);
  void buildTotalCells (char *);
  void dealWithTotalCells (char, char *, double *, double *);
  void computeTERAndSEAnalytically (int, CELL *, double *, double *);
  void computeAllMERs (int, CELL *, double *);

 public:

  CellOperations (const GetData& inputData) :
    TwoArrays (inputData), BootstrapCI (inputData)
  {
    data = &inputData;

    whichRate = data->whichRate;

    inputCISFilename = data->inputCISFilename;

    if (data->flag == 'c')
      inputCISFilename2 = data->inputCISFilename2;

    totalProb = tpSE = 0.0;
  }

  ~CellOperations ()
  {
    free (cells);
  }

};

double CellOperations::computeMERate (int WhichOne,
				    double dG, double dA, double dg, double da)
{
  double (* rateFunctions []) (double, double) =
    {
      computeWeightedRate,
      computeAverageRate
    };

  return (rateFunctions [WhichOne] (dg / dG, da / dA));
}

double CellOperations::computeWeightedRate (double fn, double fp)
{
  if ((fn + fp) != 0.0)
    return ((fn * fn + fp * fp) / (fn + fp));
  else
    return (0.0);
}

double CellOperations::computeAverageRate (double fn, double fp)
{
  return ((fn + fp) / 2.0);
}

int CellOperations::countNumCells (char *inputCISFilename)
{
  char buffer [STRINGSIZE];
  int num;
  FILE *fp;

  if ((fp = fopen (inputCISFilename, "r")) == NULL)
    {
      printf ("There is no input file for CIS.\n");
      exit (1);
    }

  num = 0;
  while (fgets (buffer, STRINGSIZE, fp) != NULL)
    {
      num++;
    }

  rewind (fp);
  fclose (fp);

  return (num - 1);
}

void CellOperations::buildCellArrays (int numCells, char *inputCISFilename,
				      CELL *cells)
{
  const char *ct0 = " \n\t=_:,";
  char * tmp, buffer [STRINGSIZE];
  int count, num;
  FILE *fp;

  if ((fp = fopen (inputCISFilename, "r")) == NULL)
    {
      printf ("There is something wrong with opening the input file.\n");
      exit (1);
    }

  if (fgets (buffer, STRINGSIZE, fp) == NULL)
    {
      printf ("There is something wrong with reading the first line of \
the input CIS file.\n");
      exit (1);
    }

  Util::bzeroOneArray (numCells, cells);

  count = 0;

  while (fgets (buffer, STRINGSIZE, fp) != NULL)
    {
      tmp = NULL;

      if ((tmp = (char *) strtok (buffer, ct0)) == NULL)
	{
	  printf ("One entry in the data file is illegal.\n");
	  exit (1);
	}

      nG (cells, count) = atoi (tmp);

      num = 0;
      while ((tmp = (char *) strtok (NULL, ct0)) != NULL)
	{
	  if (num == 0)
	    {
	      nA (cells, count) = atoi (tmp);
	    }
	  else if (num == 1)
	    {
	      na (cells, count) = atoi (tmp);
	    }
	  else if (num == 2)
	    {
	      ng (cells, count) = atoi (tmp);

	      if ((nG (cells, count) - ng (cells, count)) != 
		  (nA (cells, count) - na (cells, count)))
		{
		  printf ("There is something wrong with the %i-th data.\n",
			  count + 1);
		  exit (1);
		}

	      break;
	    }

	  num++;
	}

      Util::allocateOneArray (nG (cells, count), &gt (cells, count));
      assignValuesToGT (nG (cells, count), ng (cells, count),
			gt (cells, count));

      Util::allocateOneArray (nA (cells, count), &alg (cells, count));
      assignValuesToAlg (nA (cells, count), na (cells, count),
			 alg (cells, count));

      count++;
    }

  fclose (fp);
}

int CellOperations::determineMaxArrayLength (int numCells, CELL *cells)
{
  int i, maxLength;

  maxLength = 0;

  for (i = 0; i < numCells; i++)
    {
      if (nG (cells, i) > maxLength)
	maxLength = nG (cells, i);

      if (nA (cells, i) > maxLength)
	maxLength = nA (cells, i);
    }

  return (maxLength);
}

void CellOperations::computeTotalProbAndSE (char flag,
					    int numCells, CELL *cells,
					    double * TP, double * SE)
{
  int i;
  double denominator, weight, merate, totalprob, standError, varSum;

  denominator = totalprob = varSum = 0.0;

  for (i = 0; i < numCells; i++)
    denominator += (double) nG (cells, i);

  for (i = 0; i < numCells; i++)
    {
      weight = (double) nG (cells, i) / denominator;

      merate = computeMERate (whichRate,
			      (double) nG (cells, i), (double) nA (cells, i),
			      (double) ng (cells, i), (double) na (cells, i));

      totalprob += weight * merate;

      if (flag == '1')
	{
	  standError = bootstrapSEOfMERate ('2', whichRate,
					    nG (cells, i), nA (cells, i),
					    ng (cells, i), na (cells, i),
					    gt (cells, i), alg (cells, i),
					    computeMERate);

	  varSum += weight * weight * standError * standError;
	}
    }  

  *TP = totalprob;

  *SE = sqrt (varSum);
}

void CellOperations::buildTotalCells (char *inputCISFilename)
{
  numCells = countNumCells (inputCISFilename);

  Util::allocateOneArray (numCells, &cells);

  buildCellArrays (numCells, inputCISFilename, cells);
}

void CellOperations::dealWithTotalCells (char flag, char *inputCISFilename,
					 double *tp, double *se)
{
  buildTotalCells (inputCISFilename);

  maxArrayLength = determineMaxArrayLength (numCells, cells);

  initializeArraysBS (maxArrayLength, bsReplicationNum);

  if (flag == '1')
    computeTotalProbAndSE ('1', numCells, cells, tp, se);
}

void CellOperations::computeTERAndSEAnalytically (int numCells, CELL *cells,
						  double *ter, double *se)
{
  int i;
  double dG, dA, dg, da, denominator, weight, merate, totalprob, maxDiff,
         fn, fp, vg, va, varOne, varSum;

  denominator = totalprob = varSum = 0.0;

  for (i = 0; i < numCells; i++)
    denominator += (double) nG (cells, i);

  for (i = 0; i < numCells; i++)
    {
      dG = (double) nG (cells, i);
      dA = (double) nA (cells, i);
      dg = (double) ng (cells, i);
      da = (double) na (cells, i);

      if ((int) (dG - dg) != (int) (dA - da))
	{
	  printf ("1. The numbers are wrong.\n");
	  exit (1);
	}

      weight = dG / denominator;

      merate = computeMERate (whichRate, dG, dA, dg, da);

      totalprob += weight * merate;

      fn = dg / dG;
      fp = da / dA;

      vg = fn * (1.0 - fn) / dG;
      va = fp * (1.0 - fp) / dA;

      /*
      varOne = 0.5 * 0.5 * vg + 0.5 * 0.5 * va +
	2.0 * 0.5 * 0.5 * sqrt (vg) * sqrt (va) * 1.0;
      */
      varOne = 0.5 * sqrt (vg) + 0.5 * sqrt (va);
      varOne *= varOne;

      varSum += weight * weight * varOne;
    }

  *ter = totalprob;

  *se = sqrt (varSum);
}

void CellOperations::computeAllMERs (int numCells, CELL *cells, double *MERs)
{
  int i;

  for (i = 0; i < numCells; i++)
    {
      MERs [i] = computeMERate (whichRate,
				(double) nG (cells, i), (double) nA (cells, i),
				(double) ng (cells, i), (double) na (cells, i));
    }
}

// ========== class CIS ==========

class CIS: public CellOperations {

  const GetData * data;

  char CISFlag;
  int whichRate, nG, nA, ng, na, *GT, *Alg, iterationNum;
  double seBootstrapCIS;

  void statisticMERate ();
  void printCISMERResults (double, double, double, double, double, double,
			   double, double);

  void statisticTotalProb ();
  void printCISTPSEResults (double, double, double, double);

 public:

  CIS (const GetData& inputData) : CellOperations (inputData)
  {
    data = &inputData;

    CISFlag = data->CISFlag;
    whichRate = data->whichRate;
    if ((CISFlag == '1') || (CISFlag == 'n'))
      iterationNum = data->iterationNum;

    if (CISFlag == '1')
      {
	nG = data->nG;
	nA = data->nA;
	ng = data->ng;
	na = data->na;
      }
  }

  ~CIS ()
  {
    if (CISFlag == '1')
      {
	free (GT);
	free (Alg);
      }
  }

  void analysis_Bootstrap ();
  void analysis_AnalyticalMethod ();
  void compareMERsOfTwoAlgs ();

};

void CIS::analysis_Bootstrap ()
{
  if (CISFlag == '1')
    {
      statisticMERate ();
    }
  else if (CISFlag == 'n')
    {
      statisticTotalProb ();
    }
}

void CIS::statisticMERate ()
{
  double meRate, standError;

  createGtAndAlgArrays (nG, nA, ng, na, &GT, &Alg);

  meRate = computeMERate (whichRate,
			  (double) nG, (double) nA, (double) ng, (double) na);

  initializeArraysBS (max2 (nG, nA), bsReplicationNum);

  obtainRNGSEED ();

  standError = bootstrapSEOfMERate ('1', whichRate, nG, nA, ng, na, GT, Alg,
				    computeMERate);

  computeBCaCI (whichRate, nG, nA, ng, na, GT, Alg, meRate, computeMERate);

  printCISMERResults (meRate, standError, BCaCILo, BCaCIUp, CILo, CIUp,
		   meRate - ZSCORE * standError, meRate + ZSCORE * standError);
}

void CIS::printCISMERResults (double meRate, double seBootstrapCIS,
			      double BCaCILo, double BCaCIUp,
			      double CILo, double CIUp,
			      double ciLo, double ciUp)
{
  FILE *fp;

  if ((fp = fopen ("CISMERate", "w")) == NULL)
    {
      printf ("The output file of CISMER cannot be opened.\n");
      exit (1);
    }

  fprintf (fp, "Misclassification error rate is = %0.18f\n", meRate);
  fprintf (fp, "SE of MERate is = %0.18f\n", seBootstrapCIS);
  fprintf (fp, "CI_2 of MERate is (BCa)        = (%0.18f, %0.18f)\n",
	   BCaCILo, BCaCIUp);
  fprintf (fp, "CI_2 of MERate is (percentile) = (%0.18f, %0.18f)\n",
	   CILo, CIUp);
  fprintf (fp, "CI_2 of MERate is (normal)     = (%0.18f, %0.18f)\n",
	   ciLo, ciUp);
}

void CIS::statisticTotalProb ()
{
  int i, j;
  double denominator, weight, merate, totalprob, standError, varSum, *bsTPSE,
         CILo, dummyMedium, CIUp;
  FILE *fp, *fpR;

  obtainRNGSEED ();

  if (iterationNum == 1)
    {
      dealWithTotalCells ('1', inputCISFilename, &totalProb, &tpSE);

      printCISTPSEResults (totalProb, tpSE,
			 totalProb - ZSCORE * tpSE, totalProb + ZSCORE * tpSE);
    }
  else
    {
      Util::allocateOneArray (iterationNum, &bsTPSE);
      Util::bzeroOneArray (iterationNum, bsTPSE);

      buildTotalCells (inputCISFilename);

      maxArrayLength = determineMaxArrayLength (numCells, cells);

      initializeArraysBS (maxArrayLength, bsReplicationNum);

      denominator = totalprob = 0.0;

      for (i = 0; i < numCells; i++)
	denominator += (double) nG (cells, i);

      for (i = 0; i < numCells; i++)
	{
	  weight = (double) nG (cells, i) / denominator;

	  merate = computeMERate (whichRate,
			       (double) nG (cells, i), (double) nA (cells, i),
			       (double) ng (cells, i), (double) na (cells, i));

	  totalprob += weight * merate;
	}

      for (i = 0; i < iterationNum; i++)
	{
	  varSum = 0.0;

	  for (j = 0; j < numCells; j++)
	    {
	      weight = (double) nG (cells, j) / denominator;

	      standError = bootstrapSEOfMERate ('2', whichRate,
						nG (cells, j), nA (cells, j),
						ng (cells, j), na (cells, j),
						gt (cells, j), alg (cells, j),
						computeMERate);

	      varSum += weight * weight * standError * standError;
	    }

	  bsTPSE [i] = sqrt (varSum);
	}

      if (((fp = fopen ("bootstrapReps_SEOfTER", "w")) == NULL) ||
	  ((fpR = fopen ("__REPLICATIONS_TEMP_FILE_1__", "w")) == NULL))
	{
	  printf ("There is something wrong with opening a file.\n");
	  exit (1);
	}

      for (i = 0; i < iterationNum; i++)
	{
	  fprintf (fp, "%0.18f\n", bsTPSE [i]);
	  fprintf (fpR, "%0.18f\n", bsTPSE [i]);
	}

      fclose (fp);
      fclose (fpR);

      printf ("SE:        %0.18f\n",
	      computeUnbiasedStDev (iterationNum, bsTPSE));

      runRPercentile ("__REPLICATIONS_TEMP_FILE_1__",
		      bsSigLevel / 2.0, 1.0 - bsSigLevel / 2.0,
		      &CILo, &dummyMedium, &CIUp);

      system ("rm -f __REPLICATIONS_TEMP_FILE_1__");

      printf ("CI:       (%0.18f, %0.18f)\n", CILo, CIUp);

      qsort (bsTPSE, iterationNum, sizeof (double), doubleCompAscend);

      printf ("Min & Max: %0.18f, %0.18f\n",
	      bsTPSE [0], bsTPSE [iterationNum - 1]);
    }
}

void CIS::printCISTPSEResults (double totalProb, double tpSE,
			       double ciLo, double ciUp)
{
  FILE *fp;

  if ((fp = fopen ("CISTotalProb", "w")) == NULL)
    {
      printf ("The output file of CISTP cannot be opened.\n");
      exit (1);
    }

  fprintf (fp, "Total probability is = %0.18f\n", totalProb);
  fprintf (fp, "SE of TP is = %0.18f\n", tpSE);
  fprintf (fp, "CI_2 of TP is = (%0.18f, %0.18f)\n", ciLo, ciUp);

  fclose (fp);
}

void CIS::analysis_AnalyticalMethod ()
{
  buildTotalCells (inputCISFilename);

  computeTERAndSEAnalytically (numCells, cells, &totalProb, &tpSE);

  printf ("Total probability is = %0.18f\n", totalProb);
  printf ("SE of TP is = %0.18f\n", tpSE);
  printf ("CI_2 of TP is = (%0.18f, %0.18f)\n",
	  totalProb - ZSCORE * tpSE, totalProb + ZSCORE * tpSE);
}

void CIS::compareMERsOfTwoAlgs ()
{
  int numCells1, numCells2;
  double *MERs1, *MERs2;

  buildTotalCells (inputCISFilename);
  numCells1 = numCells;

  Util::allocateOneArray (numCells, &MERs1);
  computeAllMERs (numCells, cells, MERs1);

  buildTotalCells (inputCISFilename2);
  numCells2 = numCells;

  if (numCells1 != numCells2)
    {
      printf ("The numbers of cells in the two data files are not equal.\n");

      exit (1);
    }

  Util::allocateOneArray (numCells, &MERs2);
  computeAllMERs (numCells, cells, MERs2);

  compareTwoArrays (numCells, MERs1, MERs2);
}

// ========== main ==========

int main (int argc, char ** argv)
{
  GetData data (argc, argv);

  Util util (data);

  if (data.flag == 'b')
    {
      CIS cis (data);
      cis.analysis_Bootstrap ();
    }
  else if (data.flag == 'a')
    {
      CIS cis (data);
      cis.analysis_AnalyticalMethod ();
    }
  else if (data.flag == 'c')
    {
      CIS cis (data);
      cis.compareMERsOfTwoAlgs (); 
   }

  return 0;
}

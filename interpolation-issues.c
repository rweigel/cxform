/* To run, execute 'make interpoation-issues'*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cxform.h"

#define NUM_IGRF_YEARS_DEFINED 25

int main()
{
	int retVal, year, month, day, hour, minute, second;
	long es;
	double jd;
	double fracYear, fracYearIndex;

  /* https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt
     g/h n m 1900.0 1905.0 1910.0 1915.0 1920.0 1925.0 1930.0 1935.0 1940.0 1945.0 1950.0 1955.0 1960.0 1965.0 1970.0 1975.0 1980.0 1985.0 1990.0 1995.0   2000.0    2005.0    2010.0    2015.0   2020.0 2020-25
     g  1  0 -31543 -31464 -31354 -31212 -31060 -30926 -30805 -30715 -30654 -30594 -30554 -30500 -30421 -30334 -30220 -30100 -29992 -29873 -29775 -29692 -29619.4 -29554.63 -29496.57 -29441.46 -29404.8     5.7 
  */

	static double g01[NUM_IGRF_YEARS_DEFINED+1] =
		{-31543,    -31464,  -31354,   -31212,   -31060,    -30926,     -30805,     -30715,
		 -30654,    -30594,  -30554,   -30500,   -30421,    -30334,     -30220,     -30100,
		 -29992,    -29873,  -29775,   -29692,   -29619.4,  -29554.63,  -29496.57,  -29441.46,
	     -29404.8,  -29404.8};

	year = 2025;
	month = 12;
	day = 31;
	hour = 0;
	minute = 0;
	second = 0;

	// date2es() converts a standard Gregorian date and UT (YYYY/MM/DD HH:MM:SS) to 
	// ephemeris seconds past J2000, as required by CXFORM.
	es = date2es(year, month, day, hour, minute, second);

	// Note on fracYearIndex:
	//   157788000 = 5*365.25*86400 (365.25 is length of Julian year)
	//   -3155803200 = es on 1899-12-31 at 00:00:00
	// The fracYearIndex calculation is not exact because it assumes that each 5-year interval
	// has an average year length of 365.25 days. (In the interval 1900-1904, the average is 365.2;
	// 1900 is not a leap year and 1904 is, so 365.2 = 366 + 4*365). This approximation will
	// result in the interpolation being made for the wrong time, but the maximum error will be a fraction
	// of a day (and interpolation is done over a 5-year period). This error is probably less than
	// the uncertainty of the IRGF coefficient [https://www.ngdc.noaa.gov/IAGA/vmod/igrfhw.html].
	// However, it will lead to difficulty in comparing the resutls of coordinate transform libraries.
	//
	// There is an additional complication. The IGRF coefficients apply to 
	// 5-year periods starting in the year 1900. So it would seem that the first coefficent
	// should correspond to the midpoint of the interval [1900, 1905].
	// Here, the first coefficient applies to 1899-12-31 at 00:00:00. This will also lead to
	// a small error, probably less than the uncertainty of the IRGF coefficient.
	fracYearIndex = (es+3155803200.0)/157788000.0;
	//fracYearIndex = ((double)es+3155673600.0)/(5.0*365.2*86400.0);
	fracYear = fmod(fracYearIndex, 1.0);

	unsigned long Ne = -1 + (sizeof g01)/sizeof g01[0];

	jd = gregorian_calendar_to_jd(year, month, day, hour, minute, second);
	printf("Time: %.4d/%.2d/%.2d %.2d:%.2d:%.2d  (JD %f, ES: %ld)\n\n", 
		year, month, day, hour, minute, second, jd, es);
	fflush(stdout);

	printf("# of IGRF epochs               = %lu\n", Ne);
	printf("Last allowed date              = %lu-12-31\n", 1900+(Ne*5)-1);
	printf("fracYear (fracEpoch)           = %f\n", fracYear);
	printf("fracYearIndex (fracEpochIndex) = %f\n", fracYearIndex);
	printf("g01 at start of epoch          = %.2f\n", g01[(int)floor(fracYearIndex)]);
	printf("g01 at end of epoch            = %.2f\n", g01[(int)ceil(fracYearIndex)]);

	  if (fracYearIndex > NUM_IGRF_YEARS_DEFINED)
	  {
		fprintf(stderr, "ERROR: Specified year is greater than IGRF implementation.  Exiting.");
		exit(EXIT_FAILURE);
	  }

	
	return 0; 
}

static char const Ident[] = 
     "@(#) $Id: cxform-manual.c,v 1.21 1999/10/01 18:17:25 esm Exp $ ";
/*
** cxform-manual.c  --  manually coded functions for coordinate transform
**
** This file is part of Ed's "cxform" Coordinate Transformation package.
** It contains the hand-coded functions for converting between some
** coordinate systems.  
**
** All code in here is derived from Mike Hapgood <M.Hapgood@rl.ac.uk>'s
** excellent introduction to coordinate transformations, at
**
**	http://sspg1.bnsc.rl.ac.uk/Share/Coordinates/ct_home.htm
**
** Converted to C by Ed Santiago <esm@lanl.gov>, who takes all blame
** for any bugs you may find.
*/

#include <stdio.h>
#include <stddef.h>
#include <math.h>

#include "cxform.h"

#define X 0
#define Y 1
#define Z 2

#define DEGREES_TO_RADIANS      (M_PI / 180.0)
#define sind(x)                 sin((x) * DEGREES_TO_RADIANS)
#define cosd(x)                 cos((x) * DEGREES_TO_RADIANS)
#define atan2d(y,z)		(atan2(y,z) / DEGREES_TO_RADIANS)

#define	SECONDS_PER_CENTURY	(86400 * 365.25 * 100)


/* for debugging */
#define	DUMP_MAT	{ int i,j; for (i=0;i<3;i++) { for (j=0;j<3;j++) printf("%15lf ", mat[i][j]); printf("\n"); }}


/*****************************************************************************\
|*                                                                           *|
|*                          MATRIX OPERATIONS                                *|
|*                                                                           *|
|*  This section provides functions for defining & dealing with matrices.    *|
|*                                                                           *|
\*****************************************************************************/

/******************\
|* mat_times_mat  *|  multiplies two 3x3 matrices together.
\******************/
void mat_times_mat( Mat m1,  Mat m2, Mat m_out)
{
  int i,j;
  Mat m_tmp;

  /*
  ** Do the multiplication, but do so into a _temporary_ holding spot,
  ** in case we have been invoked with m_out == (m1 or m2).
  */
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      m_tmp[i][j] = m1[i][0]*m2[0][j] + m1[i][1]*m2[1][j] + m1[i][2]*m2[2][j];
    }
  }

  /*
  ** Now write into the caller's space.
  */
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      m_out[i][j] = m_tmp[i][j];
    }
  }
}


/******************\
|* mat_times_vec  *|  multiplies a 3x3 matrix by a 3-D vector
\******************/
void mat_times_vec(Mat m1, Vec v1, Vec v_out)
{
  int i;
  Vec v_tmp;

  for (i=0; i<3; i++) {
    v_tmp[i] = m1[i][0]*v1[0] + m1[i][1]*v1[1] + m1[i][2]*v1[2];
  }

  for (i=0; i<3; i++) {
    v_out[i] = v_tmp[i];
  }
}


/******************\
|* mat_transpose  *|  returns the transpose of a 3x3 matrix
\******************/
void mat_transpose(Mat m_in, Mat m_out)
{
  int i,j;
  Mat m_tmp;

  /*
  ** Do the transposition, but do so into a _temporary_ holding spot,
  ** in case we have been invoked with m_out == m_in
  */
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      m_tmp[i][j] = m_in[j][i];
    }
  }

  /*
  ** Now write into the caller's space.
  */
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      m_out[i][j] = m_tmp[i][j];
    }
  }
}



/*******************\
|* hapgood_matrix  *|  defines a rotation matrix for a given angle & axis
\*******************|
 *
 * Rotation matrices are a special case.  They can be defined by two 
 * parameters: an axis about which to rotate (X, Y, Z) and an angle.
 * Given those two, we can fill in all nine elements of a 3x3 matrix.
 *
 * See http://sspg1.bnsc.rl.ac.uk/Share/Coordinates/matrix.htm
 */
void hapgood_matrix(const double theta, int axis, Mat mat)
{
  int i,j;

  /* 1.calculate sin(zeta) and cos(zeta),  */
  double sin_theta = sind(theta);
  double cos_theta = cosd(theta);

  /* compute the indices for the other two axes (e.g., "X,Z" for Y) */
  int t1 = (axis+1) % 3;
  int t2 = (axis+2) % 3;
  if (t1 > t2) {
    int tmp;
    tmp = t1;
    t1  = t2;
    t2  = tmp;
  }
    

  /*
  ** 4.set the remaining off-diagonal terms to zero.
  */
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      mat[i][j] = 0.0;

  /*
  ** 2.determine the matrix diagonal:
  **   1.put 1 in the Nth term, where N=1 if the rotation axis is X, etc
  **   2.put cos(zeta) in the other two terms
  */
  mat[axis][axis] = 1.0;
  mat[t1][t1]     = cos_theta;
  mat[t2][t2]     = cos_theta;

  /*
  ** 3.locate the two off-diagonal terms in the same columns and rows as 
  **   the cos(zeta) terms - put sin(zeta) in the term above the diagonal
  **   and -sin(zeta) in the term below,
  */
  mat[t1][t2]     =  sin_theta;
  mat[t2][t1]     = -sin_theta;
}


/*****************************************************************************\
|*                                                                           *|
|*                          ELEMENTAL FUNCTIONS                              *|
|*                                                                           *|
|*  This section provides functions required by the actual transformations   *|
|*                                                                           *|
\*****************************************************************************/

/*
** H
**
** time, in hours, since preceding UT midnight
*/
double
H(const double et)
{
  double jd    = (et / 86400.0) - 0.5;	/* Convert seconds to days */
  double dfrac = jd - (long)jd;
  double hh    = dfrac * 24.0;

  if (hh < 0.0)
    hh += 24.0;

  /*  printf("hh= %15.10lf\n", hh);*/
  return hh;
}


/*
** T0
**
** Julian Centuries from a certain time to 1 Jan 2000 12:00
*/
double
T0(const double et)
{
  double jd    = (et / 86400.0) - 0.5;	/* Convert seconds to days */

  jd = (long)jd - 0.5;
  /* printf("jd= %15.10lf\n", jd); */
  return jd / 36525.0;
}


/*
** MJD
**
** Modified Julian Date
*/
double
MJD(const double et)
{
  double jd = (et / 86400.0) - 0.5;

  return jd + 51545.0;
}



/*
** lambda0
**
** The Sun's ecliptic longitude (lambdaO) can be calculated using the 
** series of formulae:
**
**     M = 357.528 + 35999.050T0 + 0.04107H degrees 
**     Lambda = 280.460 + 36000.772T0 + 0.04107H degrees 
**     lambdaO = Lambda + (1.915 - 0.0048T0) sinM + 0.020 sin2M
**
** where T0 is the time in Julian centuries from 12:00 UT on 1 January 2000
** to the midnight Universal Time (UT) preceding the time of interest and 
** H is the time in hours since that preceding UT midnight. Formulae 
** derived from the Almanac for Computers. In the intermediate formulae, 
** M is the Sun's mean anomaly and Lambda its mean longitude.
*/
double
lambda0(const double et)
{
  double M, lambda;

  M      = 357.528 + 35999.050 * T0(et) + 0.04107 * H(et);
  lambda = 280.460 + 36000.772 * T0(et) + 0.04107 * H(et);

  return lambda + (1.915 - 0.0048 * T0(et)) * sind(M) + 0.020 * sind(2 * M);
}


/*
** The obliquity of the ecliptic (epsilon) can be calculated using the formula:
**
**	epsilon = 23.439 - 0.013T0 degrees
**
** where T0 is the time in Julian centuries from 12:00 UT on 1 January 2000
** to the midnight (UT) preceding the time of interest. Formula derived from
** the Almanac for Computers.
*/
double
epsilon(const double et)
{
  return 23.439 - 0.013 * T0(et);
}


/*
** Latitude and longitude of Earth's magnetic pole
*/
double mag_lat(double et)
{
  return  78.8 + 4.283E-2 * (MJD(et) - 46066) / 365.25;
}

double mag_lon(double et)
{
  return 289.1 - 1.413E-2 * (MJD(et) - 46066) / 365.25;
}


/*****************************************************************************\
|*                                                                           *|
|*                    TRANSFORMATION MATRICES                                *|
|*                                                                           *|
|* Hapgood defines all his transformations in terms of various matrices T1,  *|
|* T2, ... Tn.  Since they're used in various places throughout the paper,   *|
|* it makes sense to define them as subroutines.                             *|
|*                                                                           *|
|* Comments preceding each function are extracted verbatim from Hapgood.     *|
|*                                                                           *|
\*****************************************************************************/

/*
** The GEI2000 to GEI transformation is given by the matrix 
**
**	P = <-zA,Z>*<thetaA,Y>*<-zetaA,Z>
**
** where the rotation angles zA, thetaA and zetaA are the precession 
** angles. This transformation is a precession correction as described
** by Hapgood (1995).
**
**     zA     = 0.64062 T0 + 0.00030 T0^2 degrees 
**     thetaA = 0.55675 T0 - 0.00012 T0^2 degrees 
**     zetaA  = 0.64062 T0 + 0.00008 T0^2 degrees
*/
void
mat_P(const double et, Mat mat)
{
  Mat mat_tmp;
  double t0 = T0(et);

  hapgood_matrix(0.64062 * t0  +  0.00030 * t0*t0, Z, mat);

  hapgood_matrix(0.55675 * t0  -  0.00012 * t0*t0, Y, mat_tmp);
  mat_times_mat(mat, mat_tmp, mat);

  hapgood_matrix(0.64062 * t0  +  0.00008 * t0*t0, Z, mat_tmp);
  mat_times_mat(mat, mat_tmp, mat);
}


/*
** The GEI to GEO transformation is given by the matrix
**
**	T1 = <theta,Z>
**
** where the rotation angle theta is the Greenwich mean sidereal time. This
** transformation is a rotation in the plane of the Earth's equator from 
** the First Point of Aries to the Greenwich meridian. 
**---
** The Greenwich sidereal time (theta) can be calculated using the formula:
**
**	theta = 100.461 + 360000.770T0 + 15.04107H degrees
**
** where T0 is the time in Julian centuries from 12:00 UT on 1 January 2000
** to the midnight Universal Time (UT) preceding the time of interest and 
** H is the time in hours since that preceding UT midnight. Formula derived
** from the Almanac for Computers.
*/
void
mat_T1(const double et, Mat mat)
{
  double theta = 100.461 + 36000.770 * T0(et) + 15.04107 * H(et);

  theta = fmod(theta, 360.0);
  if (theta < 0.0)
    theta += 360.0;

/*printf("T0= %15.10lf, H= %15.10lf, theta= %15lf\n", T0(et), H(et), theta);*/

  hapgood_matrix(theta, Z, mat);
}


/*
** The GEI to GSE transformation is given by the matrix 
**
**	T2 = <lambdaO,Z>*<epsilon,X>
**
** where the rotation angle lambdaO is the Sun's ecliptic longitude and 
** the angle epsilon is the obliquity of the ecliptic. This transformation 
** is a rotation from the Earth's equator to the plane of the ecliptic
** followed by a rotation in the plane of the ecliptic from the First Point 
** of Aries to the Earth-Sun direction. 
*/
void
mat_T2(const double et, Mat mat)
{
  Mat mat_tmp;

  hapgood_matrix(lambda0(et), Z, mat);

  hapgood_matrix(epsilon(et), X, mat_tmp);
  mat_times_mat(mat, mat_tmp, mat);
}



/*
** vec_Qe
**
** don't ask.
*/
void
vec_Qe(double et, Vec Qe)
{
  double lat = mag_lat(et);
  double lon = mag_lon(et);

  double cos_lat = cosd(lat);
  double sin_lat = sind(lat);
  double cos_lon = cosd(lon);
  double sin_lon = sind(lon);

  Mat mat_tmp, mat;

  Vec Qg;

  Qg[0] = cos_lat * cos_lon;
  Qg[1] = cos_lat * sin_lon;
  Qg[2] = sin_lat;

  /* printf("lat=%lf  lon=%lf\n", 90.0 - lat, lon);*/

  mat_T2(et, mat);
  mat_T1(et, mat_tmp);
  mat_transpose(mat_tmp, mat_tmp);
  mat_times_mat(mat, mat_tmp, mat);
  mat_times_vec(mat, Qg, Qe);
}


/*
** The GSE to GSM transformation is given by the matrix
**
**	T3 = <-psi,X>
**
** where the rotation angle psi is the the GSE-GSM angle. This
** transformation is a rotation in the GSE YZ plane from the GSE Z axis
** to the GSM Z axis. 
*/
void
mat_T3(const double et, Mat mat)
{
  Vec Qe;
  double psi;

  vec_Qe(et, Qe);

  psi = atan2d(Qe[1], Qe[2]);

  hapgood_matrix(-psi, X, mat);
}


/*
** The GSM to SM transformation is given by the matrix
**	T4 = <- mu,Y>
**
** where the rotation angle mu is the dipole tilt. This transformation 
** is a rotation in the GSM XZ plane from the GSM Z axis to the 
** geomagnetic dipole axis. 
*/
void
mat_T4(const double et, Mat mat)
{
  Vec Qe;
  double mu;

  vec_Qe(et, Qe);

  mu = atan2d(Qe[0], sqrt(Qe[1]*Qe[1] + Qe[2]*Qe[2]));

  hapgood_matrix(-mu, Y, mat);
}


/*
** The GEO to MAG transformation is given by the matrix 
**
**	T5 = <lat-90,Y>*<long,Z>
**
** where the rotation angle lat is the latitude and angle long is the 
** longitude of the geomagnetic pole (as defined by the axis of the 
** dipole component of the geomagnetic field). This transformation is 
** a rotation in the plane of the Earth's equator from the Greenwich 
** meridian to the meridian containing the dipole axis, followed by 
** a rotation in that meridian from the rotation axis to the dipole axis. 
*/
void
mat_T5(const double et, Mat mat)
{
  Mat mat_tmp;

  hapgood_matrix(mag_lat(et) - 90, Y, mat);
  hapgood_matrix(mag_lon(et),      Z, mat_tmp);
  mat_times_mat(mat, mat_tmp, mat);
}


/**************  Heliocentric transformations ************/


/*
** The HAE to HEE transformation is given by the matrix
**
**	S1 = <lambdaO + 180,Z>
**
** where the rotation angle lambdaO is the Sun's ecliptic longitude.  
** This transformation is a rotation in the plane of the ecliptic from 
** the First Point of Aries to the Sun-Earth direction. 
*/
void
mat_S1(const double et, Mat mat)
{
  hapgood_matrix(lambda0(et) + 180.0, Z, mat);
}


/*
** The HAE to HEEQ transformation is given by the matrix
**
** 	S2 = <theta0,Z>*<i,X>*<Omega,Z>
**
** where the rotation angle theta0 is the the longitude of the Sun's
** central meridian, i is the the inclination of the Sun's equator and
** Omega is the the ecliptic longitude of the ascending node of the Sun's
** equator. This transformation comprises a rotation in the plane of the
** ecliptic from the First Point of Aries to the ascending node of the
** solar equator, then a rotation from the plane of the ecliptic to the 
** plane of the equator and finally a rotation in the plane of the solar
** equator from the ascending node to the central meridian. 
*/
void
mat_S2(const double et, Mat mat)
{
  /* UNIMPLEMENTED */
}



/*****************************************************************************\
|*                                                                           *|
|*                    TRANSFORMATIONS BEGIN HERE                             *|
|*                                                                           *|
|*                                                                           *|
\*****************************************************************************/


/*********************\
**  simple_rotation  **  utility function used by all the "twixt" functions
**********************
**
** This is basically what all the "twixt" functions do:
**
**    1) define a rotation matrix
**    2) If doing an inverse transformation, transpose that matrix.
**    3) multiply the rotation matrix by the input vector.
**
** To save all that work in the "twixt" functions, they just call this
** function, passing us a pointer to a function that defines the matrix.
*/
int
simple_rotation(const double et, Vec v_in, Vec v_out, Direction d, void (*m)())
{
  Mat mat;

  /*
  ** Call the user-specified function to get a rotation matrix.
  */
  (m)(et, mat);

  /*
  ** To do the inverse transformation, we use the transposition of the matrix
  */
  if (d == BACK)
    mat_transpose(mat, mat);

  /*
  ** Multiply the rotation matrix by the input vector, and that's it!
  */
  mat_times_vec(mat, v_in, v_out);

  return 0;
}


/*  Bo-ring.  Just call simple_rotation() with the appropriate matrix. */
int
j2000_twixt_gei(const double et, Vec v_in, Vec v_out, Direction direction)
{
  return simple_rotation(et, v_in, v_out, direction, &mat_P);
}

int
gei_twixt_geo(const double et, Vec v_in, Vec v_out, Direction direction)
{
  return simple_rotation(et, v_in, v_out, direction, &mat_T1);
}

int
geo_twixt_mag(const double et, Vec v_in, Vec v_out, Direction direction)
{
  return simple_rotation(et, v_in, v_out, direction, &mat_T5);
}

int
gei_twixt_gse(const double et, Vec v_in, Vec v_out, Direction direction)
{
  return simple_rotation(et, v_in, v_out, direction, &mat_T2);
}

int
gse_twixt_gsm(const double et, Vec v_in, Vec v_out, Direction direction)
{
  return simple_rotation(et, v_in, v_out, direction, &mat_T3);
}

int
gsm_twixt_sm(const double et, Vec v_in, Vec v_out, Direction direction)
{
  return simple_rotation(et, v_in, v_out, direction, &mat_T4);
}



/*****************************************************************************\
|*              HELIOCENTRIC COORDINATE SYSTEMS                              *|
\*****************************************************************************/

/*
** Hapgood defines a transformation between GSE and HEE in his 1992
** paper (section 6), but this part isn't online.  
**
** The gist of it is, we rotate 180 degrees about Z, and then translate
** along X.
*/
int
gse_twixt_hee(const double et, Vec v_in, Vec v_out, Direction direction)
{
  Mat mat;

  hapgood_matrix(180, Z, mat);

  /*
  ** Note that there's no transposition here if the direction is "back";
  ** the operation works identically in both directions.
  */
  mat_times_vec(mat, v_in, v_out);

  /* Translate the X axis about the earth-sun distance */
  v_out[0] += (double)1.5e8;

  return 0;

  /*
  ** But we also need to add "R", a constant vector defined by
  **
  **      R = [ Rsun, 0, 0 ]
  **
  ** where 
  **
  **             r0 (1 - e^2)
  **    Rsun =   ------------
  **             1 + e cos(v)
  **
  **   r0 = 1.495985E8 km        	mean distance of the Sun from Earth.
  **
  **    e = 0.016709 - 0.0000418T0	eccentricity of the Sun's apparent
  **					orbit around the Earth.
  **
  **    w = 282.94 + 1.72 T0		longitude of perigee of that orbit
  **
  **    v = lambda0 - w			(see lambda0 above)
  */
}


int
hae_twixt_hee(const double et, Vec v_in, Vec v_out, Direction direction)
{
  return simple_rotation(et, v_in, v_out, direction, &mat_S1);
}

//
// genutahsky.c
//
// Create a Radiance description of a clear sky for any time and place
//
// Portions Copyright (c) Mark J. Stock., 2009, 2023
//
// Portions written by Mark J. Stock, markjstock@gmail.com
// Remainder written by Greg Ward
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#ifdef LIBNOVA
#include <libnova/julian_day.h>
#include <libnova/transform.h>
#include <libnova/refraction.h>
#include <libnova/apparent_position.h>
#include <libnova/solar.h>
#include <libnova/lunar.h>
#include <libnova/venus.h>
#include <libnova/jupiter.h>
#include <libnova/mars.h>
#else
#include "astronomy.h"
#endif

#include "timestruct.h"

#define STAR_THRESH -0.04	// sun z position above which no stars

// stuff from gensky.c
char *progname;
#define  PI		3.14159265358979323846
#define  DEGTORAD	0.0174532925
#define  SUNEFFICACY		208.		/* illuminant B (solar dir.) */
#undef  toupper
#define  toupper(c)     ((c) & ~0x20)   /* ASCII trick to convert case */

// M^-1 for Adobe RGB from http://www.brucelindbloom.com/Eqn_RGB_XYZ_Matrix.html
static float mi[3][3] = {2.041369, -0.969266, 0.0134474, -0.5649464, 1.8760108, -0.1183897, -0.3446944, 0.041556, 1.0154096};

// time zone characters from gensky.c
struct {
        char    zname[8];       /* time zone name (all caps) */
        float   zmer;           /* standard meridian */
} tzone[] = {
        {"YST", 135}, {"YDT", 120},
        {"PST", 120}, {"PDT", 105},
        {"MST", 105}, {"MDT", 90},
        {"CST", 90}, {"CDT", 75},
        {"EST", 75}, {"EDT", 60},
        {"AST", 60}, {"ADT", 45},
        {"NST", 52.5}, {"NDT", 37.5},
        {"GMT", 0}, {"BST", -15},
        {"CET", -15}, {"CEST", -30},
        {"EET", -30}, {"EEST", -45},
        {"AST", -45}, {"ADT", -60},
        {"GST", -60}, {"GDT", -75},
        {"IST", -82.5}, {"IDT", -97.5},
        {"JST", -135}, {"NDT", -150},
        {"NZST", -180}, {"NZDT", -195},
        {"", 0}
};


// elevation is in degrees above horizon
void getSunColor (const float elevation, const double turb, float* sunColor) {

  float f,x,y,xx,yy,zz,r,g,b;
  float luminance = 1.0;
  float altToPure,temp;

  // old way (from sunlight.c)

  // convert elevation angle to 0..1 scaled chromaticity blending
  f = 1. - sqrt(fabs(sin(elevation*DEGTORAD)));
  f = 0.;
    
  // convert elevation angle to xy chromacity
  x = 0.33 + f*0.06;
  y = 0.34 + f*0.05;
    
  // hyper-color chromaticity
  //x = 0.33 + f*0.14;
  //y = 0.34 + f*0.13;
  

  // new way (using color temperature)

  // define color temp as function of altitude
  // altitude in degrees above which color is pure 6500K (guess)
  altToPure = (turb-1.)*10.;
  temp = 6500.;

  if (fabs(elevation) < altToPure) {

    // reset the temp
    temp = 5000. + 1500.*sin(0.5*PI*elevation/altToPure);

    // debug
    //temp = 5000.;
    //fprintf(stdout,"# sun color temp %g K\n",temp);

    // now convert color temp to chromaticity
    temp = 1.e+3/temp;
    // this formula valid for 4000 < temp < 7000
    // from http://en.wikipedia.org/wiki/Standard_illuminant
    x = 0.244063 + 0.09911*temp + 2.9678*temp*temp - 4.607*temp*temp*temp;
    y = -3.*x*x + 2.87*x - 0.275;
    //fprintf(stdout,"# sun chromaticity %g %g\n",x,y);

    // convert chromacity to XYZ
    yy = luminance;
    xx = yy*x/y;
    zz = yy*(1.-x-y)/y;
  
    // convert XYZ to RGB
    //r = xx*mi[0][0] + yy*mi[0][1] + zz*mi[0][2];
    //g = xx*mi[1][0] + yy*mi[1][1] + zz*mi[1][2];
    //b = xx*mi[2][0] + yy*mi[2][1] + zz*mi[2][2];
    r = xx*mi[0][0] + yy*mi[1][0] + zz*mi[2][0];
    g = xx*mi[0][1] + yy*mi[1][1] + zz*mi[2][1];
    b = xx*mi[0][2] + yy*mi[1][2] + zz*mi[2][2];

  // or, the sun is a pure D65 illuminant
  } else {

    r = luminance;
    g = luminance;
    b = luminance;

  }

  // base color
  sunColor[0] = r;
  sunColor[1] = g;
  sunColor[2] = b;
    
  // finally, scale by gamma (not in this code)
  //sunColor[0] = exp(log(r)/2.2);
  //sunColor[1] = exp(log(g)/2.2);
  //sunColor[2] = exp(log(b)/2.2);

  //fprintf(stdout,"%.2e %.2e %.2e\n",sunColor[0],sunColor[1],sunColor[2]);
}


int writeSun (Time* time, const double inloc[3], const double turb, float *sunPos) {

  double alti, azim;	// apparent altitude, azimuth from observer on earth
  double discSize;		// solar disc size
  double appar_mag;		// apparent magnitude
  double maxAngularSegmentSize = 1.0;	// maximum size of source disc
					//  <0.53 means subsample the sun

  fprintf(stdout,"\n# Sun brightness and position from Earth's surface\n");

  // get sun position
#ifdef LIBNOVA
  struct ln_lnlat_posn obs;	// observer
  struct ln_equ_posn equ;	// equatorial sun position
  struct ln_hrz_posn hrz;	// horiz alt/az
  obs.lat = inloc[0];
  obs.lng = inloc[1];
  ln_get_solar_equ_coords (time->jd, &equ);
  ln_get_hrz_from_equ (&equ, &obs, time->jd, &hrz);
  // 360 deg azimuth is due South, 270 is due East
  // so correct it to 0=North
  alti = hrz.alt;
  azim = hrz.az+180.0;
  if (azim > 360.0) azim -= 360.0;

  // sun disc size (always about 0.53 deg)
  discSize = 2.*ln_get_solar_sdiam(time->jd)/3600.0;

  // get brightness (at ground?)
  appar_mag = −26.74;
#else

  astro_observer_t observer = Astronomy_MakeObserver(inloc[0], inloc[1], inloc[2]);
  astro_equatorial_t equ_ofdate = Astronomy_Equator(BODY_SUN, &(time->atime), observer, EQUATOR_OF_DATE, ABERRATION);
  // MUST adjust sun position due to refraction in atmosphere
  astro_horizon_t hor = Astronomy_Horizon(&(time->atime), observer, equ_ofdate.ra, equ_ofdate.dec, REFRACTION_NORMAL);
  // azimuth from Astronomy_Horizon is CW from 0=North
  alti = hor.altitude;
  azim = hor.azimuth;

  // get solar disc size in deg
  discSize = 2.0*RAD2DEG*atan(SUN_RADIUS_KM / (KM_PER_AU*equ_ofdate.dist));

  // get brightness (at ground?)
  astro_illum_t illum = Astronomy_Illumination(BODY_SUN, time->atime);
  appar_mag = illum.mag;
#endif

  fprintf(stdout,"# solar altitude %7.3f deg, azimuth %7.3f deg, size %.4f deg\n",alti,azim,discSize);
  fprintf(stdout,"# magnitude %7.4f\n",appar_mag);

  // position in vector format
  sunPos[0] = sin(azim*DEGTORAD)*cos(alti*DEGTORAD);
  sunPos[1] = cos(azim*DEGTORAD)*cos(alti*DEGTORAD);
  sunPos[2] = sin(alti*DEGTORAD);

  // sun is too low, do not draw
  if (alti < -20.0) return(false);

  // NEW

  // if sun is larger, we get more energy
  const double disc_size_multiple = pow(discSize/0.536, 2);

  // irradiance (E_e) is in W/m^2
  // at Earth distance, average solar irradiance is 1361 W/m^2 (over all spectra)
  // https://en.wikipedia.org/wiki/Solar_irradiance
  const double irrad_toa = 1361. * disc_size_multiple;

  // radiance (L_e) is in W/(sr m^2) and is what goes into the source material
  // using average solar disc size of 0.536 degrees, that's 6.87344079e-05 sr
  // https://en.wikipedia.org/wiki/Sun
  // making radiance at Earth distance 1361/6.87344079e-05 = 1.9801e+7 W/(sr m^2)
  const double rad_toa = irrad_toa / 6.87344079e-05;

  // but this 1361W is over the whole spectrum - how much are in R,G,B bands? less than 1/3rd!
  // found numbers indicating 40% of power is in visible bands

  // and assuming an equal distribution across rgb bands, that's 2640133 each, top-of-atmosphere
  const double radperband_toa = 0.4 * rad_toa / 3.;

  // illuminance (E_v) is in lux, or lumens/m^2
  // https://en.wikipedia.org/wiki/Illuminance#Relation_to_luminance
  // at peak sensitivity, 1W = 683 lumens
  // but luminous efficacy of sunlight at 5800K is 93 lm/W (or 98?)
  // https://en.wikipedia.org/wiki/Luminous_efficacy#Lighting_efficiency
  // note that 5800K truncated to the visible spectrum is 251 lm/W
  // higher because 1W in the visible spectrum is brighter than 1 W spread across a blackbody spectrum
  // making illuminance 1361 * 98 = 133378

  // another calculation for illuminance (E_v) is from "apparent magnitude"
  // https://en.wikipedia.org/wiki/Illuminance#Astronomy
  //const double illuminance = pow(10.0, 0.4*(-14.18-appar_mag));
  // Sun has an apparent magnitude of −26.74 (also given by the library)
  // making E_v = 105682 lm/m^2

  // why the difference? because the apparent magnitude is measured at sea level

  // the sky attenuates the sun's radiation
  // at Earth's surface, 1361 W/m^2 is attenuated to 1050 at most, and another 70 from the sky
  // note that 1050/1361 = 0.77149, while 105682/133378 = 0.79235
  // let's call it 22% for each normal optical depth
  // then we scale radiance by pow(0.78, 1/cos(a)), where a=0 for directly overhead

  // how to scale with altitude?
  // as a first guess, we're only scaling the integrated density of the air above us
  // conveniently, that's measured by pressure, the ratio vs. sea level can be
  // approximated using the scale height of 8400m as p/p_0 = exp(-h/8400)
  // https://en.wikipedia.org/wiki/Barometric_formula
  const double altitude_factor = 1.0 - 0.22*exp(-inloc[2] / 8400.);

  // and scale by sun elevation angle (increased optical depth)
  // but avoid division by zero and allow horizontal sun to have some brightness
  const double angle_multiple = 1.0 / sin((0.05*90. + 0.95*fmax(0.0,alti))*DEGTORAD);

  // giving the per-channel radiance at altitude given sun height as
  const double radperband_atloc = radperband_toa * pow(altitude_factor, angle_multiple);

  // so for an overhead sun on a clear day at sea level, that's 2059304
  const double lum = radperband_atloc;

  //fprintf(stderr,"# sun stuff %g %g %g %g\n", radperband_toa, altitude_factor, angle_multiple, radperband_atloc);

  // ALTERNATE
  // irradiance is in W/m^2, and needs a luminous efficacy
  //const double irradiance = pow(10.0, 0.4*(-14.18-appar_mag)) / 98.;
  // note that the magnitudes from Astronomy are not corrected for atmosphere! the are top-of-atmosphere

  // OLD
  // sun brightness (from gensky.c, color.h) at ground level
  // 1.5e+9 is 103100 / sun size in sr (from above)
  // but why is sun efficacy twice what we expect
  // and why is this not split into color bands?
  //double lum = 1.5e9/SUNEFFICACY * (1.147 - .147/(sunPos[2]>.16?sunPos[2]:.16));

  // give option of breaking sun up into many smaller suns
  // this will allow smoother penumbras at the cost of extra computation
  if (discSize > maxAngularSegmentSize) {
    // segment the disc into a large number of smaller discs,
    // this allows smoother and more accurate penumbras to be created
  }

  // someday we will handle eclipses

  // get the sun color
  float sunColor[3];
  getSunColor(alti, turb, sunColor);

  // write the Radiance sun description
  fprintf(stdout,"void light solar\n");
  fprintf(stdout,"0\n0\n3 %g %g %g\n",lum*sunColor[0],lum*sunColor[1],lum*sunColor[2]);

  // complete the Radiance sun object
  fprintf(stdout,"solar source sun\n");
  fprintf(stdout,"0\n0\n4 %g %g %g %.4f\n",sunPos[0],sunPos[1],sunPos[2],discSize);

  return(true);
}


int writeSky (float turb, float* sunPos) {

  fprintf(stdout,"\n# Sky color, luminance from Utah sky model\n");

  // dump the colorfunc
  fprintf(stdout,"void colorfunc skyfunc\n");
  fprintf(stdout,"4 skyr skyg skyb utah.cal\n0\n");
  fprintf(stdout,"4 %g %g %g %g\n",turb,sunPos[0],sunPos[1],sunPos[2]);

  // write these later
  if (false) {
    // the glow source
    fprintf(stdout,"skyfunc glow skyglow\n");
    fprintf(stdout,"0\n0\n4 1. 1. 1. 0\n");

    // and apply them both to the domes
    fprintf(stdout,"skyglow source skydome\n");
    fprintf(stdout,"0\n0\n4 0 0 1 180\n");
    fprintf(stdout,"skyglow source grounddome\n");
    fprintf(stdout,"0\n0\n4 0 0 -1 180\n");
  }

  return(true);
}


int writeMoon (Time* time, const double inloc[3]) {

  double alti, azim;	// apparent altitude, azimuth from observer on earth
  double discSize;		// solar disc size from libnova
  double adjAlt;		// altitude adjustment
  double lum;			// luminance, best guess
  double appar_mag;		// apparent magnitude
  double phase,discFrac,limb;
  float lunPos[3],lunC[3];

  // set lunar color (slightly brownish)
  lunC[0] = 1.05;
  lunC[1] = 1.00;
  lunC[2] = 0.85;
  // should really modify this near the horizon!

  // get moon position
#ifdef LIBNOVA
  struct ln_lnlat_posn obs; // observer
  struct ln_equ_posn equ;	// equatorial moon position
  struct ln_hrz_posn hrz;	// horiz alt/az

  obs.lat = inloc[0];
  obs.lng = inloc[1];

  ln_get_lunar_equ_coords (time->jd, &equ);
  ln_get_hrz_from_equ (&equ, &obs, time->jd, &hrz);
  // 360 deg azimuth is due South, 270 is due East
  // so correct it to 0=North
  alti = hrz.alt;
  azim = hrz.az+180.0;
  if (azim > 360.0) azim -= 360.0;

  // get disc size
  discSize = 2.*ln_get_lunar_sdiam(time->jd)/3600.0;

  // get phase and fraction illuminated
  appar_mag = 0.0;
  phase = ln_get_lunar_phase(time->jd);
  discFrac = ln_get_lunar_disk(time->jd);

  // get altitude adjustment due to refraction (altitude, p in millibars, temp in C)
  adjAlt = ln_get_refraction_adj (hrz.alt,1010.,10.);

  // get phase details -- later
  // limb = ln_get_lunar_bright_limb(time->jd);

#else

  astro_observer_t observer = Astronomy_MakeObserver(inloc[0], inloc[1], inloc[2]);
  astro_equatorial_t equ_ofdate = Astronomy_Equator(BODY_MOON, &(time->atime), observer, EQUATOR_OF_DATE, ABERRATION);
  // MUST adjust sun position due to refraction in atmosphere
  astro_horizon_t hor = Astronomy_Horizon(&(time->atime), observer, equ_ofdate.ra, equ_ofdate.dec, REFRACTION_NORMAL);
  // azimuth from Astronomy_Horizon is CW from 0=North
  alti = hor.altitude;
  azim = hor.azimuth;

  // get lunar disc size in deg
  discSize = 2.0*RAD2DEG*atan(MOON_EQUATORIAL_RADIUS_KM / (KM_PER_AU*equ_ofdate.dist));

  // get phase and fraction illuminated
  astro_illum_t illum = Astronomy_Illumination(BODY_MOON, time->atime);
  phase = illum.phase_angle;
  appar_mag = illum.mag;
  discFrac = illum.phase_fraction;
#endif

  // position in vector format
  lunPos[0] = sin(azim*DEGTORAD)*cos(alti*DEGTORAD);
  lunPos[1] = cos(azim*DEGTORAD)*cos(alti*DEGTORAD);
  lunPos[2] = sin(alti*DEGTORAD);

  fprintf(stdout,"\n# Lunar altitude %7.3f deg, azimuth %7.3f deg, size %6.3f deg\n",alti,azim,discSize);
  fprintf(stdout,"# magnitude %7.4f, phase %7.3f, disc illum fraction %7.3f\n",appar_mag,phase,discFrac,discSize);

  // if too low, do not draw
  if (alti < -10.0) return(false);

  // phase 0/360 is full, 180 is new

  // illuminance is in lux, or lumens/m^2
  const double illuminance = pow(10.0, 0.4*(-14.18-appar_mag));
  // but light sources use radiance, which is W/sr/m^2

  // luminance is 1/449000 of full bright sun
  // and scaled by fraction visible
  lum = 15.6 * pow(discFrac,2);

  // moon is a source, like the sun
  fprintf(stdout,"void light lunar\n");
  fprintf(stdout,"0\n0\n3 %g %g %g\n",lum*lunC[0],lum*lunC[1],lum*lunC[2]);

  // complete the moon
  fprintf(stdout,"lunar source moon\n");
  fprintf(stdout,"0\n0\n4 %g %g %g %.3f\n",lunPos[0],lunPos[1],lunPos[2],discSize);

  return(true);
}


void writePlanets (Time* time, const double inloc[3]) {

  float pos[3],col[3],mars_col[3],jup_col[3];
  double alti, azim;	// apparent altitude, azimuth from observer on earth
  double lum;			// luminance, best guess
  double discFrac,discSize;

  // white
  col[0] = 1.0;
  col[1] = 1.0;
  col[2] = 1.0;

  // martian red
  mars_col[0] = 1.2;
  mars_col[1] = 0.8;
  mars_col[2] = 0.4;

  // jupiter brown
  jup_col[0] = 1.1;
  jup_col[1] = 0.8;
  jup_col[2] = 0.5;

  // set up data structures
#ifdef LIBNOVA
  struct ln_lnlat_posn obs; // observer
  struct ln_equ_posn equ;   // equatorial moon position
  struct ln_hrz_posn hrz;   // horiz alt/az
  obs.lat = inloc[0];
  obs.lng = inloc[1];
#else
  astro_observer_t observer = Astronomy_MakeObserver(inloc[0], inloc[1], inloc[2]);
  astro_equatorial_t equ_ofdate;
  astro_horizon_t hor;
  astro_illum_t illum;
#endif

  // Venus
#ifdef LIBNOVA
  ln_get_venus_equ_coords (time->jd, &equ);
  ln_get_hrz_from_equ (&equ, &obs, time->jd, &hrz);
  // disc size should be 10-66 arcseconds
  discSize = 2.*ln_get_venus_sdiam(time->jd)/3600.0;
  // apparent magnitude
  lum = ln_get_venus_magnitude(time->jd);
  alti = hrz.alt;
  azim = hrz.az+180.0;
  if (azim > 360.0) azim -= 360.0;
#else
  equ_ofdate = Astronomy_Equator(BODY_VENUS, &(time->atime), observer, EQUATOR_OF_DATE, ABERRATION);
  hor = Astronomy_Horizon(&(time->atime), observer, equ_ofdate.ra, equ_ofdate.dec, REFRACTION_NORMAL);
  alti = hor.altitude;
  azim = hor.azimuth;
  discSize = 2.0*RAD2DEG*atan(VENUS_RADIUS_KM / (KM_PER_AU*equ_ofdate.dist));
  illum = Astronomy_Illumination(BODY_VENUS, time->atime);
  lum = illum.mag;
  discFrac = illum.phase_fraction;
#endif

  pos[0] = sin(azim*DEGTORAD)*cos(alti*DEGTORAD);
  pos[1] = cos(azim*DEGTORAD)*cos(alti*DEGTORAD);
  pos[2] = sin(alti*DEGTORAD);

  fprintf(stdout,"\n# Venus at alt %7.3f deg, az %7.3f deg, magnitude %7.3f, disc %6.2f asec\n",alti,azim,lum,3600*discSize);
  if (alti > -5.0) {
    // convert to luminance
    lum = 0.000142*exp(-0.921*lum);
    fprintf(stdout,"void light venusian\n");
    fprintf(stdout,"0\n0\n3 %g %g %g\n",lum*col[0],lum*col[1],lum*col[2]);
    fprintf(stdout,"venusian source venus\n");
    fprintf(stdout,"0\n0\n4 %g %g %g %.3f\n",pos[0],pos[1],pos[2],discSize);
  }

  // Jupiter
#ifdef LIBNOVA
  ln_get_jupiter_equ_coords (time->jd, &equ);
  ln_get_hrz_from_equ (&equ, &obs, time->jd, &hrz);
  // disc size 30-49 arcseconds
  discSize = sqrt(ln_get_jupiter_equ_sdiam(time->jd)*ln_get_jupiter_pol_sdiam(time->jd));
  discSize = 2.*discSize/3600.0;
  // apparent magnitude
  lum = ln_get_jupiter_magnitude(time->jd);
  alti = hrz.alt;
  azim = hrz.az+180.0;
  if (azim > 360.0) azim -= 360.0;
#else
  equ_ofdate = Astronomy_Equator(BODY_JUPITER, &(time->atime), observer, EQUATOR_OF_DATE, ABERRATION);
  hor = Astronomy_Horizon(&(time->atime), observer, equ_ofdate.ra, equ_ofdate.dec, REFRACTION_NORMAL);
  alti = hor.altitude;
  azim = hor.azimuth;
  discSize = 2.0*RAD2DEG*atan(JUPITER_MEAN_RADIUS_KM / (KM_PER_AU*equ_ofdate.dist));
  illum = Astronomy_Illumination(BODY_JUPITER, time->atime);
  lum = illum.mag;
  discFrac = illum.phase_fraction;
#endif

  pos[0] = sin(azim*DEGTORAD)*cos(alti*DEGTORAD);
  pos[1] = cos(azim*DEGTORAD)*cos(alti*DEGTORAD);
  pos[2] = sin(alti*DEGTORAD);

  fprintf(stdout,"\n# Jupiter at alt %7.3f deg, az %7.3f deg, magnitude %7.3f, disc %6.2f asec\n",alti,azim,lum,3600*discSize);
  if (alti > -5.0) {
    // convert to luminance
    lum = 0.000142*exp(-0.921*lum);
    fprintf(stdout,"void light jovian\n");
    fprintf(stdout,"0\n0\n3 %g %g %g\n",lum*jup_col[0],lum*jup_col[1],lum*jup_col[2]);
    fprintf(stdout,"jovian source jupiter\n");
    fprintf(stdout,"0\n0\n4 %g %g %g %.3f\n",pos[0],pos[1],pos[2],discSize);
  }

  // Mars
#ifdef LIBNOVA
  ln_get_mars_equ_coords (time->jd, &equ);
  ln_get_hrz_from_equ (&equ, &obs, time->jd, &hrz);
  // disc size 30-49 arcseconds
  discSize = 2.*ln_get_mars_sdiam(time->jd)/3600.0;
  // apparent magnitude
  lum = ln_get_mars_magnitude(time->jd);
  alti = hrz.alt;
  azim = hrz.az+180.0;
  if (azim > 360.0) azim -= 360.0;
#else
  equ_ofdate = Astronomy_Equator(BODY_MARS, &(time->atime), observer, EQUATOR_OF_DATE, ABERRATION);
  hor = Astronomy_Horizon(&(time->atime), observer, equ_ofdate.ra, equ_ofdate.dec, REFRACTION_NORMAL);
  alti = hor.altitude;
  azim = hor.azimuth;
  discSize = 2.0*RAD2DEG*atan(MARS_EQUATORIAL_RADIUS_KM / (KM_PER_AU*equ_ofdate.dist));
  illum = Astronomy_Illumination(BODY_MARS, time->atime);
  lum = illum.mag;
  discFrac = illum.phase_fraction;
#endif

  pos[0] = sin(azim*DEGTORAD)*cos(alti*DEGTORAD);
  pos[1] = cos(azim*DEGTORAD)*cos(alti*DEGTORAD);
  pos[2] = sin(alti*DEGTORAD);

  fprintf(stdout,"\n# Mars at alt %7.3f deg, az %7.3f deg, magnitude %7.3f, disc %6.2f asec\n",alti,azim,lum,3600*discSize);

  if (alti > -5.0) {
    // convert to luminance
    lum = 0.000142*exp(-0.921*lum);
    fprintf(stdout,"void light martian\n");
    fprintf(stdout,"0\n0\n3 %g %g %g\n",lum*mars_col[0],lum*mars_col[1],lum*mars_col[2]);
    fprintf(stdout,"martian source mars\n");
    fprintf(stdout,"0\n0\n4 %g %g %g %.3f\n",pos[0],pos[1],pos[2],discSize);
  }

  // no other plants are ever brighter than Mars (really?)
}


int writeStars (Time* time, double zPos, const double inloc[3]) {

  // some positions

  float vnorth[3],veq[3];
  float sinth,costh,rz1,rx1,rz2;

  // scale the brightness of the star map
  float starBright = 3.e-4;

  // note: largest stars are .06 arcseconds in diameter
  //       typical stars are .007 arcseconds in diameter
  //       sky (360 deg) is 1.3M arcseconds wide
  //       thus, all stars should be one pixel in any practical resolution

  // are any of the stars bright enough to see?
  // is sun higher than 2 degrees below the horizon?
  if (zPos > STAR_THRESH) return(false);

  // otherwise, we see stars

#ifdef LIBNOVA
  struct ln_lnlat_posn obs; // observer
  obs.lat = inloc[0];
  obs.lng = inloc[1];
  struct ln_equ_posn astar;
  struct ln_equ_posn zeromotion;
  struct ln_equ_posn equ;
  struct ln_hrz_posn hrz;	// horiz alt/az
  // to figure out how to position the dome, use
  // http://libnova.sourceforge.net/group__apparent.html
  // (void) ln_get_apparent_posn (mean_position, proper_motion, time->jd, apparent_position);
  // do that once for a pole star, and once for a equatorial star, and convert!

  // set proper motion to 0,0 because it won't affect much
  zeromotion.ra = 0.0;
  zeromotion.dec = 0.0;

  // hypothetical pole star
  astar.ra = 0.0;
  astar.dec = 90.0;
  // find vector to pole
  ln_get_apparent_posn (&astar,&zeromotion,time->jd,&equ);
  //fprintf(stderr,"north star ra %g  dec %g\n",equ.ra,equ.dec);
  ln_get_hrz_from_equ (&equ, &obs, time->jd, &hrz);

  vnorth[0] = -sin(hrz.az*DEGTORAD)*cos(hrz.alt*DEGTORAD);
  vnorth[1] = -cos(hrz.az*DEGTORAD)*cos(hrz.alt*DEGTORAD);
  vnorth[2] = sin(hrz.alt*DEGTORAD);
  //fprintf(stderr,"north star az %g  alt %g\n",hrz.az,hrz.alt);
  //fprintf(stderr,"north star x,y,z %g %g %g\n",vnorth[0],vnorth[1],vnorth[2]);

  // find second rotation (azimuth is always ~180)
  rx1 = hrz.alt - 90.;
  // find third rotation
  rz2 = hrz.az - 180.;

  // hypothetical equatorial star
  astar.ra = 0.0;
  astar.dec = 0.0;
  // find vector to zero-meridian on equator
  ln_get_apparent_posn (&astar,&zeromotion,time->jd,&equ);
  //fprintf(stderr,"equatorial star ra %g  dec %g\n",equ.ra,equ.dec);
  ln_get_hrz_from_equ (&equ, &obs, time->jd, &hrz);

  veq[0] = -sin(hrz.az*DEGTORAD)*cos(hrz.alt*DEGTORAD);
  veq[1] = -cos(hrz.az*DEGTORAD)*cos(hrz.alt*DEGTORAD);
  veq[2] = sin(hrz.alt*DEGTORAD);
  //fprintf(stderr,"equatorial star az %g  alt %g\n",hrz.az,hrz.alt);
  //fprintf(stderr,"equatorial star x,y,z %g %g %g\n",veq[0],veq[1],veq[2]);

#else
  // set observer and test stars
  astro_observer_t observer = Astronomy_MakeObserver(inloc[0], inloc[1], inloc[2]);

  // hypothetical equatorial star
  // args are body, ra (0..24), dec (-90..90), dist (LY)
  Astronomy_DefineStar(BODY_STAR1, 0.0, 0.0, 1.0e+6);
  // hypothetical pole star
  Astronomy_DefineStar(BODY_STAR2, 0.0, 90.0, 1.0e+6);

  // vector to zero-meridian on equator
  astro_equatorial_t equ_ofdate = Astronomy_Equator(BODY_STAR1, &(time->atime), observer, EQUATOR_OF_DATE, ABERRATION);
  astro_horizon_t hor = Astronomy_Horizon(&(time->atime), observer, equ_ofdate.ra, equ_ofdate.dec, REFRACTION_NONE);
  veq[0] = sin(hor.azimuth*DEGTORAD)*cos(hor.altitude*DEGTORAD);
  veq[1] = cos(hor.azimuth*DEGTORAD)*cos(hor.altitude*DEGTORAD);
  veq[2] = sin(hor.altitude*DEGTORAD);
  //fprintf(stderr,"equatorial star az %g  alt %g\n",hor.azimuth,hor.altitude);
  //fprintf(stderr,"equatorial star x,y,z %g %g %g\n",veq[0],veq[1],veq[2]);

  // vector to pole
  equ_ofdate = Astronomy_Equator(BODY_STAR2, &(time->atime), observer, EQUATOR_OF_DATE, ABERRATION);
  hor = Astronomy_Horizon(&(time->atime), observer, equ_ofdate.ra, equ_ofdate.dec, REFRACTION_NONE);
  vnorth[0] = sin(hor.azimuth*DEGTORAD)*cos(hor.altitude*DEGTORAD);
  vnorth[1] = cos(hor.azimuth*DEGTORAD)*cos(hor.altitude*DEGTORAD);
  vnorth[2] = sin(hor.altitude*DEGTORAD);
  //fprintf(stderr,"north star az %g  alt %g\n",hor.azimuth,hor.altitude);
  //fprintf(stderr,"north star x,y,z %g %g %g\n",vnorth[0],vnorth[1],vnorth[2]);

  // find second rotation (azimuth is always ~180)
  rx1 = hor.altitude - 90.;
  // find third rotation
  rz2 = hor.azimuth;// - 180.;

#endif

  // first z-rotation
  // NOTE: this first rotation needs a fixed offset, but I don't know
  // what to use. Maybe render a view at midnight and then run Stallarium
  // to see how much to shift?
  sinth = veq[2]/sin(rx1*DEGTORAD);
  if (sinth < -1.) sinth = -1.;
  if (sinth > 1.) sinth = 1.;
  //fprintf(stderr,"sinth %g\n",sinth);
  rz1 = asinf(sinth);
  //fprintf(stderr,"theta %g or %g\n",rz1*57.2957787,180-rz1*57.2957787);
  // now, do we use this, or pi-this ?
  costh = (veq[0] + sin(rz2*DEGTORAD)*cos(rx1*DEGTORAD)*sinth) / cos(rz2*DEGTORAD);
  //fprintf(stderr,"costh %g\n",costh);
  if (costh < 0.) rz1 = 180. - (rz1/DEGTORAD);
  else rz1 = rz1/DEGTORAD;

  // include the proper file with the proper rotations
  fprintf(stdout,"\n# Stars from Tycho-2 star catalog\n");
  fprintf(stdout,"!xform -rz %g -rx %g -rz %g stardome.rad\n",rz1,rx1,rz2);

  // old way: create the geometry calls right here
  if (false) {
    fprintf(stdout,"\n# Stars from Tycho-2 star catalog\n");
 
    // create the imagemap
    fprintf(stdout,"void colorpict starmapcolor\n");
    fprintf(stdout,"7 noneg noneg noneg TychoSkymapII.t5_08192x04096.hdr sphere.cal inf_u inf_v\n");
    fprintf(stdout,"0\n1 0.5\n");

    // then the glow source
    fprintf(stdout,"starmapcolor glow starmapglow\n");
    fprintf(stdout,"0\n0\n4 %g %g %g -1\n",starBright,starBright,starBright);

    // and the source
    fprintf(stdout,"starmapglow source starmap\n");
    fprintf(stdout,"0\n0\n4 0 0 1 360\n");

    // what about a low-resolution light map?
    fprintf(stdout,"void colorpict starcolor\n");
    fprintf(stdout,"7 noneg noneg noneg TychoSkymapII.t5_00080x00040.hdr sphere.cal inf_u inf_v\n");
    fprintf(stdout,"0\n1 0.5\n");

    // then the glow source
    // should this be an "illum"?
    //fprintf(stdout,"\nstarcolor glow starglow\n");
    //fprintf(stdout,"0\n0\n4 %g %g %g 0\n",starBright,starBright,starBright);
    fprintf(stdout,"starcolor illum starglow\n");
    fprintf(stdout,"1 void\n0\n3 %g %g %g\n",starBright,starBright,starBright);

    // and the source
    fprintf(stdout,"starglow source starlight\n");
    fprintf(stdout,"0\n0\n4 0 0 1 360\n");

    // now, we need to mix these! No, we don't.
    //fprintf(stdout,"\nvoid mixfunc stars\n");
    //fprintf(stdout,"0\n0\n4 0 0 1 360\n");
  }

  return(true);
}


// write the mixfunc to merge the sky and the stars
void writeSkyStarMix (double zPos) {

  if (zPos > STAR_THRESH) {
    // no stars, write sky as usual
    fprintf(stdout,"\n# Applying sky color map\n");

    // the glow source
    fprintf(stdout,"skyfunc glow skyglow\n");
    fprintf(stdout,"0\n0\n4 1. 1. 1. 0\n");

    // and apply it to the domes
    fprintf(stdout,"skyglow source skydome\n");
    fprintf(stdout,"0\n0\n4 0 0 1 180\n");
    fprintf(stdout,"skyglow source grounddome\n");
    fprintf(stdout,"0\n0\n4 0 0 -1 180\n");

  } else {
    // mix the stars and sky
    fprintf(stdout,"\n# Mixing sky and star color maps\n");

    // mixing the sky colors
    fprintf(stdout,"void mixfunc mixedcolor\n");
    fprintf(stdout,"4 skyfunc starmapcolor half half.cal\n0\n0\n");

    // glow needs 2x multiplier because of the mix
    fprintf(stdout,"mixedcolor glow mixedglow\n");
    fprintf(stdout,"0\n0\n4 2. 2. 2. 0\n");

    // and apply the mixed one to the upper dome
    fprintf(stdout,"mixedglow source skydome\n");
    fprintf(stdout,"0\n0\n4 0 0 1 190\n");
    fprintf(stdout,"mixedglow source grounddome\n");
    fprintf(stdout,"0\n0\n4 0 0 -1 170\n");
  }
}


void writeClouds (void) {
  // if geometry is available, place some clouds
}


void writeHaze (int isSun, int isSky, int isMoon, int isStars) {

  int numDirect;

  // are there enough direct light sources to color the haze?
  //numDirect = isSun + isSky + isMoon + isStars;
  numDirect = isSun + isSky + isMoon;

  // write the haze description
  if (numDirect > 0) {
    fprintf(stdout,"\n# ground haze, scale is 1unit = 33m\n");
    fprintf(stdout,"void mist hazemat\n");
    fprintf(stdout,"%d",numDirect);
    if (isSun) fprintf(stdout," sun");
    if (isSky) fprintf(stdout," skydome");
    if (isMoon) fprintf(stdout," moon");
    // stars now included in skydome!
    //if (isStars) fprintf(stdout," starmap");
    fprintf(stdout,"\n",isSun);
    // medium density
    //fprintf(stdout,"0\n7 1.2e-3 1.5e-3 1.7e-3 1 1 1 0.4\n");
    // lower density
    fprintf(stdout,"0\n7 0.6e-3 0.7e-3 0.8e-3 1 1 1 0.4\n");
  }
}


// convert hour string
int cvthour( char  *hs, double *hour, int *tsolar, double *s_meridian) {

  register char  *cp = hs;
  register int    i, j;

  if ( (*tsolar = *cp == '+') ) cp++;              /* solar time? */
  while (isdigit(*cp)) cp++;
  if (*cp == ':')
    *hour = atoi(hs) + atoi(++cp)/60.0;
  else {
    *hour = atof(hs);
    if (*cp == '.') cp++;
  }
  while (isdigit(*cp)) cp++;
  if (!*cp)
    return(0);
  if (*tsolar || !isalpha(*cp)) {
    fprintf(stderr, "%s: bad time format: %s\n", progname, hs);
    exit(1);
  }
  i = 0;
  do {
    for (j = 0; cp[j]; j++)
      if (toupper(cp[j]) != tzone[i].zname[j])
        break;
    if (!cp[j] && !tzone[i].zname[j]) {
      *s_meridian = tzone[i].zmer;
      return(1);
    }
  } while (tzone[i++].zname[0]);

  fprintf(stderr, "%s: unknown time zone: %s\n", progname, cp);
  fprintf(stderr, "Known time zones:\n\t%s", tzone[0].zname);
  for (i = 1; tzone[i].zname[0]; i++)
       fprintf(stderr, " %s", tzone[i].zname);
  putc('\n', stderr);
  exit(1);
}


/* print command header */
void printhead( register int  ac, register char  **av) {
  putchar('#');
  while (ac--) {
    putchar(' ');
    fputs(*av++, stdout);
  }
  putchar('\n');
}


/* print usage error and quit */
void userror(char  *msg) {
  if (msg != NULL)
    fprintf(stderr, "%s: Use error - %s\n", progname, msg);
  fprintf(stderr, "Usage: %s month day hour [options]\n", progname);
  exit(1);
}


int main(int argc, char **argv) {

  char errmsg[128];

  // time stuff
  Time time;		// structure holds libnova and/or astronomy date
  int year,month,day;
  int ihour,iminute;
  double dhours,dminutes,dseconds;

  // location stuff
  double observer[3];	// lat (deg), long (deg), height (m)
  bool got_meridian = false;
  double s_meridian = 0.0;

  // astronomy stuff
  float sunPos[3];
  int tsolar;
  bool isSun = false;
  bool isSky = false;
  bool isMoon = false;
  bool isStars = false;

  // atmosphere stuff
  float turbidity = 2.45;	// gensky default, also near europe average
  double gprefl = 0.2;		// deciduous forest

  // set defaults ---------------------------------------

  // set to Boston
  observer[0] = 42.36;
  observer[1] = -71.06;	// note: east longitude (use west for cli)
  observer[2] = 10.0;	// this is in meters

  set_to_now(&time);
  //write_utc(time);

  // pull year out of UTC and set as default
  year = get_year(time);

  // parse command-line arguments -----------------------
  // 1st arg is month number
  // 2nd arg is day
  // 3rd arg is 24-hour decimal hours (can add "EST" or other time zone to string!)
  // -a latitude (North assumed)
  // -o longitude (West assumed)

  // use code from gensky.c

  progname = argv[0];
  if (argc < 4)
    userror("arg count");
  month = atoi(argv[1]);
  if (month < 1 || month > 12)
    userror("bad month");
  day = atoi(argv[2]);
  if (day < 1 || day > 31)
    userror("bad day");
  got_meridian = (bool)cvthour(argv[3], &dhours, &tsolar, &s_meridian);
  for (int i = 4; i < argc; i++)
    if (argv[i][0] == '-' || argv[i][0] == '+')
      switch (argv[i][1]) {
      case 'y':
        year = atoi(argv[++i]);
        break;
      case 't':
        turbidity = atof(argv[++i]);
        break;
      case 'g':
        gprefl = atof(argv[++i]);
        break;
      case 'a':
		// keep in degrees!
        observer[0] = atof(argv[++i]);
        break;
      case 'o':
		// note negative to match gensky behavior!
        observer[1] = -atof(argv[++i]);
        break;
      case 'e':
		// elevation in meters
        observer[2] = atof(argv[++i]);
        break;
      case 'm':
        if (got_meridian) {
          ++i;
          break;          /* time overrides */
        }
		// keep in degrees
        s_meridian = atof(argv[++i]);
        break;
      default:
        sprintf(errmsg, "unknown option: %s", argv[i]);
        userror(errmsg);
      }
    else
      userror("bad option");

  // if no meridian, assume local time?
  if (!got_meridian) s_meridian = -observer[1];
  // check for meridian far away from observer location
  if (fabs(s_meridian+observer[1]) > 45.)
          fprintf(stderr,
  "%s: warning: %.1f hours btwn. standard meridian and longitude\n",
                  progname, (-observer[1]-s_meridian)/15.);

  printhead(argc, argv);

  //fprintf(stdout,"tsolar %d\n",tsolar);
  //fprintf(stdout,"got_meridian %d\n",got_meridian);
  //fprintf(stdout,"s_meridian %g\n",s_meridian);
  //fprintf(stdout,"observer.lng %g\n",observer[1]);

  fprintf(stdout,"# observer at %g deg lat, %g deg long, %g m ASL\n", observer[0], observer[1], observer[2]);

  // convert the time -----------------------------------

  ihour = floor(dhours);
  dminutes = 60.*(dhours-(double)(ihour));
  iminute = floor(dminutes);
  dseconds = 60.*(dminutes-(double)(iminute));
  set_to_given(&time, year, month, day, ihour, iminute, dseconds, s_meridian);
  write_utc(time);

  // write output ---------------------------------------

#ifdef LIBNOVA
  fprintf(stdout,"# using libnova for astronomical data\n");
#else
  fprintf(stdout,"# using Astronomy for astronomical data\n");
#endif

  // compute and write the sun "light" and "source" descriptions
  isSun = writeSun (&time, observer, turbidity, sunPos);

  // always write the sky dome
  isSky = writeSky (turbidity, sunPos);

  // place the moon if it is above the horizon
  isMoon = writeMoon (&time, observer);

  // place the three brightest planets if they are above the horizon
  writePlanets (&time, observer);

  // if the sky is dim enough, place stars
  isStars = writeStars (&time, sunPos[2], observer);
  (void) writeSkyStarMix (sunPos[2]);

  // if there is cloud geometry available, place it
  //writeClouds ();

  // create a material for the ground haze
  writeHaze (isSun,isSky,isMoon,isStars);

  // ANSI C requires main to return int
  return 0;
}

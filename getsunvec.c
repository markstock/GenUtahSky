//
// getsunvec.c
//
// Just dump the vector to the sun, assuming +x east, +y north, +z up
// and brightness equal to top-of-atmosphere
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
#else
#include "astronomy.h"
#endif

#include "timestruct.h"

// stuff from gensky.c
char *progname;
#define  PI		3.14159265358979323846
#define  DEGTORAD	0.0174532925
#define  SUNEFFICACY		208.		/* illuminant B (solar dir.) */
#undef  toupper
#define  toupper(c)     ((c) & ~0x20)   /* ASCII trick to convert case */


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


int writeSun (Time* time, const double inloc[3], const double turb, float *sunPos) {

  double alti, azim;	// apparent altitude, azimuth from observer on earth
  double discSize;		// solar disc size
  double maxAngularSegmentSize = 1.0;	// maximum size of source disc
					//  <0.53 means subsample the sun

  fprintf(stdout,"\n# Sun brightness and position from space\n");

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

  fprintf(stdout,"# libnova: solar altitude %7.3f deg, azimuth %7.3f deg, size %.4f deg\n",alti,azim,discSize);

#else
  astro_observer_t observer = Astronomy_MakeObserver(inloc[0], inloc[1], inloc[2]);
  astro_equatorial_t equ_ofdate = Astronomy_Equator(BODY_SUN, &(time->atime), observer, EQUATOR_OF_DATE, NO_ABERRATION);
  astro_horizon_t hor = Astronomy_Horizon(&(time->atime), observer, equ_ofdate.ra, equ_ofdate.dec, REFRACTION_NONE);
  // azimuth from Astronomy_Horizon is CW from 0=North
  alti = hor.altitude;
  azim = hor.azimuth;

  // get solar disc size in deg
  discSize = 2.0*RAD2DEG*atan(SUN_RADIUS_KM / (KM_PER_AU*equ_ofdate.dist));

  fprintf(stdout,"# astronomy: solar altitude %7.3f deg, azimuth %7.3f deg, size %.4f deg\n",alti,azim,discSize);
#endif

  // position in vector format
  sunPos[0] = sin(azim*DEGTORAD)*cos(alti*DEGTORAD);
  sunPos[1] = cos(azim*DEGTORAD)*cos(alti*DEGTORAD);
  sunPos[2] = sin(alti*DEGTORAD);

  // do not adjust sun position due to refraction in atmosphere

  // sun brightness (from gensky.c, color.h)
  // this is for Earth's surface, not top-of-atmosphere
  //double lum = 1.5e9/SUNEFFICACY * (1.147 - .147/(sunPos[2]>.16?sunPos[2]:.16));
  double lum = 1.5e9/SUNEFFICACY;

  // is this top-of-atmosphere?
  lum = 1.00262e+07;

#ifdef LIBNOVA
  // use solar disc size?
  //double appar_mag = ln_get_venus_magnitude(jd);
  //lum = 0.000142*exp(-0.921*appar_mag);
#else
  // get apparent magnitude of sun
  astro_illum_t sun_illum = Astronomy_Illumination(BODY_SUN, time->atime);
  double appar_mag = sun_illum.mag;
  // see https://en.wikipedia.org/wiki/Illuminance#Astronomy
  double illuminance = pow(10, 0.4*(-14.18-appar_mag));
#endif

  // give option of breaking sun up into many smaller suns
  // this will allow smoother penumbras at the cost of extra computation
  if (discSize > maxAngularSegmentSize) {
    // segment the disc into a large number of smaller discs,
    // this allows smoother and more accurate penumbras to be created
  }

  // someday we will handle eclipses

  // write the Radiance sun description
  fprintf(stdout,"void light solar\n");
  fprintf(stdout,"0\n0\n3 %g %g %g\n",lum,lum,lum);

  // complete the Radiance sun object
  fprintf(stdout,"solar source sun\n");
  fprintf(stdout,"0\n0\n4 %g %g %g %.4f\n",sunPos[0],sunPos[1],sunPos[2],discSize);

  return(true);
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

  // atmosphere stuff
  float turbidity = 2.45;	// gensky default, also near europe average

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
      case 'a':
		// keep in degrees!
        observer[0] = atof(argv[++i]);
        break;
      case 'o':
		// note negative to match gensky behavior!
        observer[1] = -atof(argv[++i]);
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

  // write output ---------------------------------------

#ifdef LIBNOVA
  fprintf(stdout,"# using libnova for atronomical data\n");
#else
  fprintf(stdout,"# using Astronomy for atronomical data\n");
#endif

  // compute and write the sun "light" and "source" descriptions
  isSun = writeSun (&time, observer, turbidity, sunPos);

  // ANSI C requires main to return int
  return 0;
}

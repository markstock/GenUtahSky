//
// getsunvec.c
//
// Just dump the vector to the sun, assuming +x east, +y north, +z up
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
#include <libnova/julian_day.h>
#include <libnova/transform.h>
#include <libnova/refraction.h>
#include <libnova/apparent_position.h>
#include <libnova/solar.h>
#include <libnova/lunar.h>
#include <libnova/venus.h>
#include <libnova/jupiter.h>
#include <libnova/mars.h>
#define TRUE 1
#define FALSE 0
#define STAR_THRESH -0.04	// sun z position above which no stars

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



int writeSun (double jd, struct ln_lnlat_posn obs, double turb, float *sunPos) {

  struct ln_equ_posn equ;	// equatorial sun position
  struct ln_hrz_posn hrz;	// horiz alt/az
  double lum;			// luminance, from gensky.c
  double discSize;		// solar disc size from libnova
  double angularSegmentSize = 1.0;	// maximum size of source disc
					//  <0.53 means subsample the sun

  // get sun position
  ln_get_solar_equ_coords (jd, &equ);
  ln_get_hrz_from_equ (&equ, &obs, jd, &hrz);
  // 360 deg azimuth is due South, 270 is due East

  // position in vector format
  sunPos[0] = -sin(hrz.az*DEGTORAD)*cos(hrz.alt*DEGTORAD);
  sunPos[1] = -cos(hrz.az*DEGTORAD)*cos(hrz.alt*DEGTORAD);
  sunPos[2] = sin(hrz.alt*DEGTORAD);

  // do not adjust sun position due to refraction in atmosphere

  fprintf(stdout,"\n# Sun brightness and position from space\n");
  fprintf(stdout,"# solar altitude %7.3f deg, azimuth %7.3f deg\n",hrz.alt,hrz.az);

  // sun brightness (from gensky.c, color.h)
  lum = 1.5e9/SUNEFFICACY * (1.147 - .147/(sunPos[2]>.16?sunPos[2]:.16));

  // sun disc size (always about 0.53 deg)
  discSize = 2.*ln_get_solar_sdiam(jd)/3600.0;

  // give option of breaking sun up into many smaller suns
  // this will allow smoother penumbras at the cost of extra computation
  if (discSize > angularSegmentSize) {
    // segment the disc into a large number of smaller discs,
    // this allows smoother penumbras to be created
  }

  // write the sun description
  fprintf(stdout,"void light solar\n");
  fprintf(stdout,"0\n0\n3 1.00262e+07 1.00262e+07 1.00262e+07\n");

  // complete the sun
  fprintf(stdout,"solar source sun\n");
  fprintf(stdout,"0\n0\n4 %g %g %g %.3f\n",sunPos[0],sunPos[1],sunPos[2],discSize);

  return(TRUE);
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
                        //*s_meridian = tzone[i].zmer * (PI/180);
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

  int i;
  char errmsg[128];

  // time stuff
  double simJD;		// Julian day, simulated
  struct ln_date date;	// date structure
  struct ln_zonedate zdate;	// date structure, includes gmtoff (seconds east of UTC)
  int year,month,day;
  double hour;

  // location stuff
  struct ln_lnlat_posn observer;
  int got_meridian = 0;
  double s_meridian = 0.0;

  // astronomy stuff
  float sunPos[3];
  int tsolar;
  int isSun = FALSE;
  int isSky = FALSE;
  int isMoon = FALSE;
  int isStars = FALSE;

  // atmosphere stuff
  float turbidity = 2.45;	// gensky default, also near europe average
  double gprefl = 0.2;		// deciduous forest

  // set defaults ---------------------------------------

  // set to my house in Newton Center
  observer.lat = 42.36;
  observer.lng = -71.06;	// note: east longitude (use west for cli)

  // set the julian date to right now
  simJD = ln_get_julian_from_sys();

  //(void) ln_get_local_date (simJD, &date);
  //fprintf(stdout,"# %4d-%02d-%02d %2d:%02d:%02d\n",
  //        date.years,date.months,date.days,
  //        date.hours,date.minutes,(int)date.seconds);

  // create the time string for the current date
  (void) ln_get_date (simJD, &date);
  //fprintf(stdout,"# UTC now: %4d-%02d-%02d %2d:%02d:%02d\n",
  //        date.years,date.months,date.days,
  //        date.hours,date.minutes,(int)date.seconds);

  // 14400 for 4 hours behind UTC
  //(void) ln_date_to_zonedate (&date, &zdate, -14400);
  //fprintf(stdout,"# Local time now: %4d-%02d-%02d %2d:%02d:%02d\n",
  //        zdate.years,zdate.months,zdate.days,
  //        zdate.hours,zdate.minutes,(int)zdate.seconds);

  // and use those as defaults (really only to set the year)
  year = date.years;
  month = date.months;
  day = date.days;
  hour = date.hours + date.minutes/60. + date.seconds/3600.;

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
        got_meridian = cvthour(argv[3], &hour, &tsolar, &s_meridian);
        for (i = 4; i < argc; i++)
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
                                //observer.lat = atof(argv[++i]) * (PI/180);
                                observer.lat = atof(argv[++i]);
                                break;
                        case 'o':
				// note negative to match gensky behavior!
                                //observer.lng = -atof(argv[++i]) * (PI/180);
                                observer.lng = -atof(argv[++i]);
                                break;
                        case 'm':
                                if (got_meridian) {
                                        ++i;
                                        break;          /* time overrides */
                                }
				// keep in degrees
                                //s_meridian = atof(argv[++i]) * (PI/180);
                                s_meridian = atof(argv[++i]);
                                break;
                        default:
                                sprintf(errmsg, "unknown option: %s", argv[i]);
                                userror(errmsg);
                        }
                else
                        userror("bad option");

        // if no meridian, assume local time?
        if (!got_meridian) s_meridian = -observer.lng;
        // check for meridian far away from observer location
        if (fabs(s_meridian+observer.lng) > 45.)
                fprintf(stderr,
        "%s: warning: %.1f hours btwn. standard meridian and longitude\n",
                        progname, (-observer.lng-s_meridian)/15.);

        printhead(argc, argv);

        //fprintf(stdout,"tsolar %d\n",tsolar);
        //fprintf(stdout,"got_meridian %d\n",got_meridian);
        //fprintf(stdout,"s_meridian %g\n",s_meridian);
        //fprintf(stdout,"observer.lng %g\n",observer.lng);


  // convert the time -----------------------------------

  // change date to match input
  //date.years = year;
  //date.months = month;
  //date.days = day;
  //date.hours = floor(hour);
  //simJD = ln_get_julian_day(&date);
  //fprintf(stdout,"# jd %15.10e\n",simJD);

  zdate.years = year;
  zdate.months = month;
  zdate.days = day;
  zdate.hours = floor(hour);
  zdate.minutes = floor(60.*(hour-(double)(zdate.hours)));
  zdate.seconds = 3600.*(hour-(double)(zdate.hours)-(double)(zdate.minutes)/60.);
  zdate.gmtoff = 86400 - (long)(240.*s_meridian);
  fprintf(stdout,"# Using time: %4d-%02d-%02d %2d:%02d:%02d\n",
          zdate.years,zdate.months,zdate.days,
          zdate.hours,zdate.minutes,(int)zdate.seconds);
  simJD = ln_get_julian_local_date(&zdate);
  fprintf(stdout,"# Julian day: %15.10e\n",simJD);

  // write output ---------------------------------------

  // compute and write the sun "light" and "source" descriptions
  isSun = writeSun (simJD, observer, turbidity, sunPos);

  // ANSI C requires main to return int
  return 0;
}

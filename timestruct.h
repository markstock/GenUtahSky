//
// timestruct.h
//
// Some routines to allow using time in both systems
//
// Copyright (c) Mark J. Stock., 2023
//

#include <stdio.h>
#include <stdlib.h>

#ifdef LIBNOVA
#include <libnova/julian_day.h>
#else
#include "astronomy.h"
#endif

// use struct to contain one or both times
typedef struct time_structure {
#ifdef LIBNOVA
  double jd;			// used in libnova
#else
  astro_time_t atime;	// used in astronomy
#endif
} Time;

void set_to_now(Time* t) {
#ifdef LIBNOVA
  t->jd = ln_get_julian_from_sys();
#else
  t->atime = Astronomy_CurrentTime();
#endif
}

void write_utc(Time t) {
#ifdef LIBNOVA
  struct ln_date date;	// date structure
  ln_get_date (t.jd, &date);
  fprintf(stdout,"# UTC: %4d-%02d-%02d %2d:%02d:%02d\n",
          date.years,date.months,date.days,
          date.hours,date.minutes,(int)date.seconds);
#else
  astro_utc_t utc = Astronomy_UtcFromTime(t.atime);
  fprintf(stdout,"# UTC: %4d-%02d-%02d %2d:%02d:%02d\n",
          utc.year,utc.month,utc.day,
          utc.hour,utc.minute,(int)utc.second);
#endif
}

int get_year(Time t) {
  int thisyear;
#ifdef LIBNOVA
  struct ln_date date;  // date structure
  ln_get_date (t.jd, &date);
  thisyear = date.years;
#else
  astro_utc_t utc = Astronomy_UtcFromTime(t.atime);
  thisyear = utc.year;
#endif
  return thisyear;
}

void set_to_given(Time* t, const int year, const int month, const int day,
                  const int hour, const int minute, const double seconds,
                  const double meridian) {
  fprintf(stdout,"# local time: %4d-%02d-%02d %2d:%02d:%02d\n",
          year,month,day,hour,minute,(int)seconds);
#ifdef LIBNOVA
  struct ln_zonedate zdate;	// date structure, includes gmtoff (seconds east of UTC)
  zdate.years = year;
  zdate.months = month;
  zdate.days = day;
  zdate.hours = hour;
  zdate.minutes = minute;
  zdate.seconds = seconds;
  //zdate.gmtoff = 86400 - (int)(240.*meridian);
  zdate.gmtoff = -(int)(240.*meridian);
  t->jd = ln_get_julian_local_date(&zdate);
#else
  astro_utc_t utc;
  utc.year = year;
  utc.month = month;
  utc.day = day;
  utc.hour = hour;
  utc.minute = minute;
  utc.second = seconds;
  astro_time_t unadj = Astronomy_TimeFromUtc(utc);
  // now adjust it based on meridian - this can't be correct...what if it changes day?
  const double days = meridian / 360.0;
  t->atime = Astronomy_AddDays(unadj, days);
#endif
}


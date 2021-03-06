#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include "gday.h"
#include "constants.h"




/* utilities */
double round_to_value(double, double);
double day_length(int, int, double);
void   calculate_daylength(state *, int, double);
int    is_leap_year(int);
void   prog_error(const char *, const unsigned int);
bool   float_eq(double, double);
void calc_warmest_quarter_temp(control *, params *, met_arrays *, 
                                met *, state *, double);

char   *rstrip(char *);
char   *lskip(char *);
char   *find_char_or_comment(char*, char);
char   *strncpy0(char*, char*, size_t);
char   *strip_first_and_last_character(char);

#endif /* UTILITIES_H */

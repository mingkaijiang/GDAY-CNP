
#include "utilities.h"



int is_leap_year(int yr) {

    if (yr % 400 == 0 || (yr % 100 != 0 && yr % 4 == 0)) {
        return TRUE;
    } else {
        return FALSE;
    }
}



double day_length(int doy, int num_days, double latitude) {

    /*

    Daylength in hours

    Eqns come from Leuning A4, A5 and A6, pg. 1196

    Reference:
    ----------
    Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.

    Parameters:
    -----------
    doy : int
        day of year, 1=jan 1
    yr_days : int
        number of days in a year, 365 or 366
    latitude : float
        latitude [degrees]

    Returns:
    --------
    dayl : float
        daylength [hrs]

    */
    double deg2rad, latr, sindec, a, b;

    deg2rad = M_PI / 180.0;
    latr = latitude * deg2rad;
    sindec = -sin(23.5 * deg2rad) * cos(2.0 * M_PI * (doy + 10.0) / num_days);
    a = sin(latr) * sindec;
    b = cos(latr) * cos(asin(sindec));

    return 12.0 * (1.0 + (2.0 / M_PI) * asin(a / b));
}

void calculate_daylength(state *s, int num_days, double latitude) {
    /* wrapper to put the day length into an array */
    int i;
    for (i = 0; i < num_days; i++) {
        s->day_length[i] = day_length(i+1, num_days, latitude);
    }
    return;
}

void prog_error(const char *reason, const unsigned int line)
{
    fprintf(stderr, "%s, failed at line: %d\n", reason, line);
	exit(EXIT_FAILURE);

    return;
}

bool float_eq(double a, double b) {
    /*
    Are two floats approximately equal...?

    Reference:
    ----------
    D. E. Knuth. The Art of Computer Programming. Sec. 4.2.2 pp. 217-8.
    */
    return fabs(a - b) <= EPSILON * fabs(a);
}



char* rstrip(char* s)
{
    /* Strip whitespace chars off end of given string, in place. Return s. */

    char* p = s + strlen(s);
    while (p > s && isspace(*--p)) *p = '\0';
    return s;
}


char* lskip(char* s)
{
    /* Return pointer to first non-whitespace char in given string. */

    while (*s && isspace(*s)) s++;
    return (char*)s;
}

char* find_char_or_comment(char* s, char c)
{
    /*

    Return pointer to first char c or ';' comment in given string, or
    pointer to null at end of string if neither found. ';' must be
    prefixed by a whitespace character to register as a comment.

    */

    int was_whitespace = 0;
    while (*s && *s != c && !(was_whitespace && *s == ';')) {
        was_whitespace = isspace(*s);
        s++;
    }
    return (char*)s;
}


char *strncpy0(char* dest, char* src, size_t size)
{
    /* Version of strncpy that ensures dest (size bytes) is null-terminated. */

    strncpy(dest, src, size);
    dest[size - 1] = '\0';
    return dest;
}

double round_to_value(double number, double roundto) {
    return (round(number / roundto) * roundto);
}


void calc_warmest_quarter_temp(control *c, params *p, met_arrays *ma, 
                                 met *m, state *s, double yr) {
    // calculate mean temperature of the warmest quarter
    // Ref: Atkin et al. (2015) New Phytologist
    // Table 6 best model for area-based broadleaved tree leaf respiration
    //
    // Parameters:
    // ----------
    // TWQ: float
    //      mean temperature of the warmest quarter;
    // q1 - q4: float
    //      quarterly-based mean air temperature;
    double TWQ, q1, q2, q3, q4;
    double t1, t2, t3, t4;
    int doy = 0;

    if (is_leap_year(yr)) {
        for (doy = 0; doy < 366; doy ++) {
            if (doy >= 0 && doy <= 90) {
                t1 += ma->tair[doy];
            } else if (doy >= 91 && doy <= 181) {
                t2 += ma->tair[doy];
            } else if (doy >= 182 && doy <= 273) {
                t3 += ma->tair[doy];
            } else if (doy >= 274 && doy < 366) {
                t4 += ma->tair[doy];
            }
        }
    } else {
        for (doy = 0; doy < 365; doy ++) {
            if (doy >= 0 && doy <= 89) {
                t1 += ma->tair[doy];
            } else if (doy >= 90 && doy <= 180) {
                t2 += ma->tair[doy];
            } else if (doy >= 181 && doy <= 272) {
                t3 += ma->tair[doy];
            } else if (doy >= 273 && doy < 365) {
                t4 += ma->tair[doy];
            }
        }
    }

    /* compute quarterly mean temperature */
    if (is_leap_year(yr)) {
       q1 = t1 / 91.0;
    } else {
       q1 = t1 / 90.0;
    }

    q2 = t2 / 91.0;
    q3 = t3 / 92.0;
    q4 = t4 / 92.0;

    TWQ = MAX(q1, MAX(q2, MAX(q3, q4)));

//    fprintf(stderr, "t1 %f, t2 %f, q1 %f, q2 %f, TWQ %f\n", 
//                     t1, t2, q1, q2, TWQ);

    s->twq = TWQ;

    return;

}

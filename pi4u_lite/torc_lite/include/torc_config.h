/* include/torc_config.h.  Generated from torc_config.h.in by configure.  */
/* include/torc_config.h.in.  Generated from configure.ac by autoheader.  */


#ifndef _PS_CONFIG
#define _PS_CONFIG


/* explicit cache line size */
#define CACHE_LINE_SIZE 128

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## _

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <pthread.h> header file. */
#define HAVE_PTHREAD_H 1

/* Define to 1 if you have the `sched_yield' function. */
#define HAVE_SCHED_YIELD 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Name of package */
#define PACKAGE "Torc-Lite"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#define PACKAGE_NAME "Torc-Lite"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "Torc-Lite 0.1a"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "Torc-Lite"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.1a"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Runtime statistics */
#define TORC_STATS /**/

/* Version number of package */
#define VERSION "0.1a"

/* not found */
/* #undef sched_yield */

#endif

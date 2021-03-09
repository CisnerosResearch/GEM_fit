#ifndef _fortran_namespace_h_
#define _fortran_namespace_h_

#ifdef NO_SECOND_UNDERSCORE
#define US2(x) x ## _
#else
#define US2(x) x ## __
#endif

#endif//_fortran_namespace_h_

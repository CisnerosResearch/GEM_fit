
#if defined __EFC
integer LOC

#elif defined __IFC
integer LOC

#elif defined __GNUC__
intrinsic IArgC
intrinsic LOC

#elif ! defined MAKEDEP
#error System type not recognized!
#endif


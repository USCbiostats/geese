#ifndef GEESE_FLOCK_CASES
#define GEESE_FLOCK_CASES 1

#define IF_GEESE(a) if (Rf_inherits((a), "geese"))
#define IF_FLOCK(a) else if (Rf_inherits((a), "flock"))
#define IF_NEITHER() else stop("The passed object is neither a 'geese' nor a 'flock'.");

#define CHECK_GEESE(a) if (!Rf_inherits((a), "geese")) \
  stop("The passed object is not of class 'geese'");

#define CHECK_FLOCK(a) if (!Rf_inherits((a), "flock")) \
  stop("The passed object is not of class 'flock'");

#endif


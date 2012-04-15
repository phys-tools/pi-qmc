#include "Hungarian.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

Hungarian::Hungarian(int size)
    :   size(size), kindex(new int[size]), work(new int[6*size]) {
}

Hungarian::~Hungarian() {
    delete [] kindex;
    delete [] work;
}

int Hungarian::operator[](int index) const {
    return kindex[index];
}

double Hungarian::getSum() const {
    return sum;
}

void Hungarian::solve(double* matrix) {
    assndx(&MODE,matrix,&size,&size,&size,kindex,&sum,work,&size);
    for (int index = 0; index < size; ++index) {
        kindex[index] -= 1;
    }
}

const int Hungarian::MODE = 1;

// assndx.f translated by f2c, then some rewriting.
void Hungarian::assndx(const int *mode, double *a, int *n,
        int *m, int *ida, int *k, double *sum, int *iw,
        int *idw) {

    static int c__1 = 1;

    /* System generated locals */
    int a_dim1, a_offset, iw_dim1, iw_offset, i__1, i__2;
    double d__1, d__2;

    /* Local variables */
    static int i__, j, j1, icl, irl, ipp, new__, irs, jsv;
    static bool lsw;
    static int icl0, icbl, imin, imax;
    static double rmin;
    static int iflag;
    static char errtxt[80];


    /* Parameter adjustments */
    a_dim1 = *ida;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --k;
    iw_dim1 = *idw;
    iw_offset = 1 + iw_dim1;
    iw -= iw_offset;

    /* Function Body */
    if (*n < 1 || *m < 1) {
        return;
    }
    imax = (*n>*m) ? *n : *m;
    imin = (*n<*m) ? *n : *m;
    *sum = 0.;
    if (*n <= *m) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    rmin = a[i__ + a_dim1];
	    i__2 = *m;
	    for (j = 1; j <= i__2; ++j) {
/* L2: */
/* Computing MIN */
		d__1 = rmin, d__2 = a[i__ + j * a_dim1];
		rmin = (d__1 < d__2) ? d__1 : d__2;
	    }
	    *sum += rmin;
	    i__2 = *m;
	    for (j = 1; j <= i__2; ++j) {
/* L3: */
		a[i__ + j * a_dim1] -= rmin;
	    }
/* L1: */
	}
    }
    if (*n >= *m) {
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    rmin = a[j * a_dim1 + 1];
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L5: */
/* Computing MIN */
		d__1 = rmin, d__2 = a[i__ + j * a_dim1];
		rmin = (d__1 < d__2) ? d__1 : d__2;
	    }
	    *sum += rmin;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L7: */
		a[i__ + j * a_dim1] -= rmin;
	    }
/* L4: */
	}
    }
    i__1 = imax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k[i__] = 0;
/* L8: */
	iw[i__ + iw_dim1] = 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    if (a[i__ + j * a_dim1] + iw[j + iw_dim1] == 0.) {
		k[i__] = j;
		iw[j + iw_dim1] = i__;
		goto L12;
	    }
/* L13: */
	}
L12:
	;
    }
L10:
    iflag = *n;
    irl = 0;
    icl = 0;
    irs = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iw[i__ + iw_dim1 * 5] = 0;
	if (k[i__] == 0) {
	    ++irl;
	    iw[irl + iw_dim1 * 6] = i__;
	    iw[i__ + iw_dim1 * 5] = -1;
	    --iflag;
	}
/* L11: */
    }
    if (iflag == imin) {
	if (*mode == 2) {
	    i__1 = imax;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L70: */
		k[i__] = iw[i__ + iw_dim1];
	    }
	}
	return;
    }
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
/* L14: */
	iw[j + (iw_dim1 << 2)] = 0;
    }
L30:
    i__ = iw[irs + iw_dim1 * 6];
    ++irs;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	if (a[i__ + j * a_dim1] + iw[j + (iw_dim1 << 2)] == 0.) {
	    iw[j + (iw_dim1 << 2)] = i__;
	    ++icl;
	    iw[icl + (iw_dim1 << 1)] = j;
	    new__ = iw[j + iw_dim1];
	    if (new__ == 0) {
		j1 = j;
L61:
		iw[j1 + iw_dim1] = iw[j1 + (iw_dim1 << 2)];
		i__ = iw[j1 + (iw_dim1 << 2)];
		if (k[i__] == 0) {
		    k[i__] = j1;
		    goto L10;
		}
		jsv = j1;
		j1 = k[i__];
		k[i__] = jsv;
		goto L61;
	    }
	    ++irl;
	    iw[irl + iw_dim1 * 6] = new__;
	    iw[new__ + iw_dim1 * 5] = i__;
	}
/* L31: */
    }
    if (irs <= irl) {
	goto L30;
    }
    lsw = true;
    icl0 = icl;
    icbl = 0;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	if (iw[j + (iw_dim1 << 2)] == 0) {
	    ++icbl;
	    iw[icbl + iw_dim1 * 3] = j;
	}
/* L41: */
    }
    rmin = a[iw[iw_dim1 * 6 + 1] + iw[iw_dim1 * 3 + 1] * a_dim1];
    i__1 = irl;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = icbl;
	for (j = 1; j <= i__2; ++j) {
/* L42: */
/* Computing MIN */
	    d__1 = rmin, d__2 = a[iw[i__ + iw_dim1 * 6] + iw[j + iw_dim1 * 3] 
		    * a_dim1];
	    rmin = (d__1 < d__2) ? d__1 : d__2; //min(d__1,d__2);
	}
    }
    *sum += rmin * (irl + icbl - imax);
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (iw[i__ + iw_dim1 * 5] == 0) {
	    i__1 = icl0;
	    for (ipp = 1; ipp <= i__1; ++ipp) {
/* L49: */
		a[i__ + iw[ipp + (iw_dim1 << 1)] * a_dim1] += rmin;
	    }
	    goto L44;
	}
	i__1 = icbl;
	for (ipp = 1; ipp <= i__1; ++ipp) {
	    new__ = iw[ipp + iw_dim1 * 3];
	    a[i__ + new__ * a_dim1] -= rmin;
	    if (lsw && a[i__ + new__ * a_dim1] + iw[new__ + (iw_dim1 << 2)] ==
		     0.) {
		iw[new__ + (iw_dim1 << 2)] = i__;
		if (iw[new__ + iw_dim1] == 0) {
		    j1 = new__;
		    lsw = false; //FALSE_;
		} else {
		    ++icl;
		    iw[icl + (iw_dim1 << 1)] = new__;
		    ++irl;
		    iw[irl + iw_dim1 * 6] = iw[new__ + iw_dim1];
		}
	    }
/* L45: */
	}
L44:
	;
    }
    if (lsw) {
	i__2 = icl;
	for (i__ = icl0 + 1; i__ <= i__2; ++i__) {
/* L51: */
	    iw[iw[iw[i__ + (iw_dim1 << 1)] + iw_dim1] + iw_dim1 * 5] = iw[i__ 
		    + (iw_dim1 << 1)];
	}
	goto L30;
    } else {
L62:
	iw[j1 + iw_dim1] = iw[j1 + (iw_dim1 << 2)];
	i__ = iw[j1 + (iw_dim1 << 2)];
	if (k[i__] == 0) {
	    k[i__] = j1;
	    goto L10;
	}
	jsv = j1;
	j1 = k[i__];
	k[i__] = jsv;
	goto L62;
    }
}

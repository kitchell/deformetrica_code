/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the MIT License. This file is also distributed     *
*    under the terms of the Inria Non-Commercial License Agreement.                    *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef DEFORMETRICA_APPROXIMATION_H
#define DEFORMETRICA_APPROXIMATION_H


#ifdef USE_FAST_MATH
	static double (*math_exp)(double x) = fast_math::fast_exp;
#else
    static double (*math_exp)(double x) = exp;
#endif


namespace fast_math {

    /**
     * Source: http://theoval.cmp.uea.ac.uk/~gcc/publications/pdf/nc2000a.pdf
     */
static double EXP_A = (1048576/M_LN2);
static double EXP_C = (60801);

inline
double fast_exp(double x) {

    union {
        double d;
#if defined(LITTLE_ENDIAN) || defined(__LITTLE_ENDIAN)
        struct { int j, i; } n;
#else
        struct { int i, j; } n;
#endif
    } _eco;

    _eco.n.i = (int)(EXP_A*(x)) + (1072693248 - EXP_C);
    _eco.n.j = 0;

    return _eco.d;
}


}

#endif //DEFORMETRICA_APPROXIMATION_H

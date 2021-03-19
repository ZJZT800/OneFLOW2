#ifndef UTILROUTINEDEFINITIONS
#define UTILROUTINEDEFINITIONS


/* *********************************************************************************
 * @file util.cpp
 * @class ArrayUtils
 * @author Kelly Black <kjblack@gmail.com>
 * @version 0.1
 * @copyright BSD 2-Clause License
 *
 * @section LICENSE
 *
 * Copyright (c) 2014, Kelly Black
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 *
 * Class to provide a set of basic utilities that are used by numerous
 * other classes.
 *
 * This is the code file for the ArrayUtils class. It includes the
 * code associated with the methods that are used to construct and
 * delete arrays used in a variety of other classes.
 *
 *
 * @brief Code file for the basic utilities associated with managing
 * arrays.
 *
 * ********************************************************************************* */


#include<iostream>
#include <cstdlib>
#include "Utils.h"

 using namespace std;
/** ************************************************************************
 * Template for allocating a two dimensional array.
 *
 *	@param n1 Number of entries for the first dimension.
 *  @param n2 Number of entries for the second dimension.
 *  @return A pointer to the array.
 * ************************************************************************ */
template <class number>
number** ArrayUtils<number>::twotensor(int n1, int n2) {
    number** u;
    register int s, i;

    //u = new number*[n1];
    u = (number**)malloc(n1 * sizeof(number*));
    u[0] = (number*)malloc(n1 * n2 * sizeof(number));

    if (u == NULL) {
        std::cout << "Error - twotensor. Could not allocate memory." << std::endl;
        std::exit(2);
    }

    if (u[0] == NULL) {
        std::cout << "Error - twotensor. Could not allocate memory for vector." << std::endl;
        std::exit(2);
    }
    for (s = 0; s < n1; s++)
    {
        u[s] = u[0] + s * n2;
    }
    for (s = 0; s < n1; s++)
    {
        for (i = 0; i < n2; i++)
        {
            u[s][i] = 0.0;
        }
    }

    return(u);
}


/** ************************************************************************
 * Template for allocating a one dimensional array.
 *
 *	@param n1 Number of entries for the dimension.
 *  @return A pointer to the array.
 * ************************************************************************ */
template <class number>
number* ArrayUtils<number>::onetensor(int n1) {
    number* u;
    register int i;

    //u = new number[n1];
    u = (number*)malloc(n1 * sizeof(number));
    if (u == NULL) {
        std::cout << "Error - onetensor. Could not allocate memory." << std::endl;
        std::exit(2);
    }

    for (i = 0; i < n1; ++i)
        u[i] = 0.0;

    return(u);

}

template <class number>
void ArrayUtils<number>::deltwotensor(number** u) {


    if (u == NULL)
        return;
    free(u[0]);
    u[0] = NULL;
    free(u);
    u = NULL;
}

template <class number>
void ArrayUtils<number>::delonetensor(number* u) {

    if (u == NULL)
        return;
    free(u);
    u = NULL;
}

#endif

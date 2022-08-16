// qprop_v3.2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "Imag.h"
#include "Real.h"
#include "tsurff.h"
#include "isurff.h"
#include "tsurff_MPI.h"

int main()
{
    //ImaginaryProp();
    //RealProp();
    // isurfv(); //ensure tsurff-version long 2 in tsurff.param
    //tsurff();  //ensure tsurff-version long 1 in tsurff.param
    tsurff_MPI();

}





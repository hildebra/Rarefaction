#ifndef randed_h
#define randed_h


/*
created using:NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
Copyright (C) 1988-1992 by Cambridge University Press. Programs Copyright (C) 1988-1992 by Numerical Recipes Software.
*/
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI 3.141592654

typedef unsigned long type;


#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>


class RAND{
public:
	RAND(long num);
	~RAND() {}
	float rand() {return randq(); }
	float rand1();
	float randq();
	float poidev(float xm);
	double gasdev();
	type bindev(double, type ,int opt = 0);

	//helper
	float gammln(float xx);
protected:
	float bind(float pp, type n);

	long idum ;
	int iset;
	double gset;

};

#endif
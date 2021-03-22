#pragma once
#include "Utils.h"

class BasicCompute
{
public:
	BasicCompute();
	~BasicCompute();

public:
	void CSRMVs(double* A, int* IA, int* JA, double* b, double* ans, int rows);
	void apxy(double* b, double* c, int rows, double realnum);     //b = realnum * b + b * c
	void apxy(double* b, double* c, int rows);                   //b = b + b * c
	double dot(double* b, double* c, int rows);
	double norm(double* b, int rows);
	void VCs(double* b, double* c, double realnum, int rows);
	void vecPlus(double* b, double* c, int rows);
	void vecPlus(double* b, double* c, int rows, double realnum);
	void vecPlus(double* b, double* c, double* d, int rows, double realnum);
	void vecPlusb(double* b, double* c, int rows, double realnum);
	void vecCopy(double* b, double* c, int rows);    //set b = c
	void rotg(double* X, double* Y, double* c, double* s);    //givens rotation; calculate coefficient S & C
	void rot(double* X, double* Y, double* c, double* s);
};
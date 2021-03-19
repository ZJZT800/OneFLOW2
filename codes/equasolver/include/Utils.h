#ifndef UTILROUTINE
#define UTILROUTINE

template <class number>
class ArrayUtils
{

public:

	/** ************************************************************************
	 * Base constructor  for the ArrayUtils class.
	 *
	 * There is not anything to do so this is an empty method.
	 *
	 *  ************************************************************************ */
	ArrayUtils() {};

	static number** twotensor(int n1, int n2);
	static number* onetensor(int n1);

	static void deltwotensor(number** u);
	static void delonetensor(number* u);

};


#include "Utils.cpp"


#endif

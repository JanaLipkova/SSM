/* 
 *
 *  Test about using blitz array
 *  and its conversion to double
 *
 *
 */

//#include <iostream>
//#include <math.h>
//#include <stdio.h>
//#include <assert.h>
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <sys/time.h>
//#include <stdarg.h>

//#include <blitz/blitz.h>
#include <blitz/array.h>
#include <iostream>
#include <stdlib.h>


using namespace std;
using namespace blitz;

int main()
{
	Array<float,2> A(3,3), B(3,3), C(3,3);
	
	A = 1, 0, 0, 
		2, 2, 2, 
		1, 0, 0;
	
	B = 0, 0, 7, 
		0, 8, 0, 
		9, 9, 9;
	
	C = A + B;
	cout << "A = " << A << endl;
	cout << "B = " << B << endl;
	cout << "C = " << C << endl;
	
	Array<double, 1> tmp(1,1), tmp2(1,1); 
	tmp = 1.1;
	tmp2(0) = 10.;
	
	
	double out = tmp(0) + tmp2(0);
	
	cout << "tmp =" << tmp << endl;
	cout << "out =" << out << endl;
	
	return 0;
}
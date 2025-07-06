%module pycrysfml
%{
#include "*.h"
%}
%include "cpointer.i"
%pointer_class(int, intp);
%pointer_class(double, doublep);
%pointer_class(float, floatp);
%include "std_string.i"
%include "cstring.i"
%include "std_vector.i"
namespace std {
	%template(FloatVector) vector<float>;
	%template(FloatMatrix) vector< vector<float> >;
	%template(IntVector) vector<int>;
}
%include "*.h"

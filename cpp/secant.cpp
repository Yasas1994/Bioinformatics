#include <iostream>
#include <cmath>
using namespace std;
double f(double x){
	
	
	return (x*x*x)+(5*x)-3;

}
double x_3(double x1, double x2){

	return x2 - ((f(x2)*(x2-x1))/(f(x2)-f(x1)));
	
	
}
 

int main(){
	// calculates roots using secant method first derivative should be given as the input in order to find the minima of the function
	double input1_x1 = 2;
	double input2_x2 = 4;
	double out_x3 = 0;
	
	
	while (abs(f(input2_x2)) > 0) {
		
		out_x3=(x_3(input1_x1,input2_x2));

		input1_x1 = input2_x2;
		input2_x2 = out_x3;

		printf ("%.25f ",f(out_x3));
		printf ("%.25f \n", out_x3);
	        
		
	}

	//cout << f(input1_x1)<<endl;	
	printf("%.25f \n",f(5.74166));
	

	return 0;
}

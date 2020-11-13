#include <iostream>
#include <cmath>

using namespace std;

double f(float x){
	
	
	return (x*x*x)+(5*x)-3;
	


}

float gradient (float x1){
	
// approximating the derivative of the function
	float h = 0.0001; 

	return (f(x1 + h) - f(x1)) / h ;

	

}

float out_x2 (float x1){

	return x1 - (f(x1)/(gradient(x1)));
}

int main(){
//Newton Raphson root finding method first derivate of the function should be given as input
	double input_x1 = 4;

	double output_x2;
	
	cout << "f(x)               x"<<endl;
	while(abs(f(input_x1))> 0.00001 ){

	output_x2 = out_x2(input_x1);
	input_x1 = output_x2;

	 //printf ("%.25f ",f(input_x1));
	 //printf("%.25f\n",input_x1);
	
	cout << std::fixed<< f(input_x1)<<" "<< input_x1<<endl;
	}



return 0;

}

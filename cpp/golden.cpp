#include <iostream>
#include <cmath>
//THIS CODE APPROXIMATES MINIMA OF A FUNCTION USING GOLDEN SECTION METHOD
using namespace std;

double f(float x){
	
	
	return (pow(x,2))-(4*x)-(10);
	

}


int main(){
	
	float x1;
	float x2;
	float a = -5;
	float b = 10;
	float phi  = 0.618033;
	
	cout<<"x1            x2             f(x1)         f(x2)\n";
//specifying precision of the output
	cout.precision(10);
	while (abs(f(a)-f(b))>0.01){

	 x2 = a + phi*(b-a);
	 x1 = b - phi*(b-a);
	
		if(f(x1)> f(x2)){
			a = x1;

		}
		if(f(x2) > f(x1)){
			b = x2;
		}
	cout<<std::fixed<<a<<" "<<b<<" "<<f(a)<<" "<<f(b)<<endl;


	

	}
	cout << "\nbest estimate: "<< (x1+x2)/2 << " "<<f((x1+x2)/2)<<endl;




return 0;
}

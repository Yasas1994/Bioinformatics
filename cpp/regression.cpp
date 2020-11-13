#include <iostream>
#include <cmath>
using namespace std;



double h(float x, float w0, float w1){

	 

	return float (w0 + (w1*x));
	}

int main(){

float y[6] = {2.0,4.0,6.0,8.0,10.0,12.0};
float x[6] = {1.0,2.0,3.0,4.0,5.0,6.0};
float w0 = 1.0;
float w1 = 1.0;
float alpha = 0.001;

//gradient calculation
float len = float(sizeof(y))/float(4.0);

for (int j = 0; j < 20000000; j++ ){
	float gradient_w0 = 0;
	float gradient_w1 = 0;

	for (int i = 0 ;i < len;i++ ){
		float e = h(x[i],w0,w1);
		gradient_w0 += (e-y[i]);
		gradient_w1 += ((e-y[i])*x[i]);
		
		}
	//cout <<gradient_w1<<endl;
	float w1new = w1 - (alpha*gradient_w1);
	float w0new = w0 - (alpha*gradient_w0);
	w1 = w1new;
	w0 = w0new;
	

	}

cout<<w1<<" "<<w0<<endl;

return 0;
}

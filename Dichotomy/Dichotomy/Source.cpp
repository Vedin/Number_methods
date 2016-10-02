#include<iostream>
#include<cmath>

using namespace std;

double func(double x) {
	return x*x*x - 3 * x*x - 14 * x - 8.0;
}
double d_func(double x) {
	return 3 * x*x - 6 * x - 14.0;
}
double EPS = 1e-7;
double dixotom() {
	double l = -10, r = 10;
	int i = 0;
	while (l + EPS < r) {
		cout << "Iter: " << i++ << endl;
		double mid = (l + r) / 2;
		cout <<fixed<< "l: " << l << " r: " << r << " mid:" << mid<<" f(mid): "<<func(mid)<<endl;
		if (abs(func(mid)) < EPS) return mid;
		if (func(mid)*func(l)<0) r = mid;
		else l = mid;
	}
	return l;
}
const double X0 = 4;
double Newton() {
	double x0 = X0;
	int i = 0;
	while (fabs(func(x0)) > EPS) {
		cout <<fixed<< "Iter: " << i << " x" << i++ << ": " << x0 << " f(x)=" << func(x0)<<endl;
		x0 = x0 - func(x0) / d_func(x0);
	}
	return x0;
}
double lyambd(double x) {
	return -0.05;
}
double fi_func(double x) {
	return x-lyambd(x)*func(x);
}
double Iter() {
	double x0 = X0;
	double xn = X0;
	int i = 0;
	do{

		cout <<fixed<< "Iter: " << i << " x" << i++ << ": " << x0 << " f(x)=" << func(x0) << endl;
		x0 = xn;

		xn = fi_func(x0);
	} while (fabs(x0 - xn) > EPS);
	return xn;
}

int main() {
	freopen("output.txt", "w", stdout);
	
	cout << "Dixotom:" << endl;
	cout<< dixotom()<<endl;
	cout << "Newton:" << endl;
	cout<< Newton() << endl;
	cout << "Simple iteration:" << endl;
	cout<< Iter() << endl;
	return 0;
}
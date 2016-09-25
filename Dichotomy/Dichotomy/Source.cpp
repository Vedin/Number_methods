#include<iostream>
#include<cmath>

using namespace std;

double func(double x) {
	return x*x*x - 3 * x*x - 14 * x - 8.0;
}
double EPS = 1e-8;
double dixotom() {
	double l = -100, r = 100;
	while (l + EPS < r) {
		double mid = (l + r) / 2;
		if (abs(func(mid)) < EPS) return mid;
		if (func(mid)*func(l)<0) r = mid;
		else l = mid;
	}
	return l;
}
int main() {
	cout << dixotom();
	return 0;
}
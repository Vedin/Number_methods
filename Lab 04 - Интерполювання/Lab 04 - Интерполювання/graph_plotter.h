#ifndef GRAPH_PLOTTER
#define GRAPH_PLOTTER

#include "polinom.h"
#include <vector>
#include <Windows.h>
using namespace std;

#define Point pair<double, double>
#define GREEN RGB(0, 255, 0)
#define RED RGB(255, 0, 0)
#define BLUE RGB(0, 0, 255)

void drawLine(HDC &hdc, Point a, Point b, COLORREF col = RGB(255, 255, 255)) {
	
	int dx = b.first - a.first;
	int dy = b.second - a.second;

	if (abs(dx) < abs(dy)) { // Primary axis Y
		if (dy < 0) swap(a, b);
		int ax = a.first, bx = b.first;
		int ay = a.second, by = b.second;
		
		if (dx == 0) {
			for (int y = ay; y <= by; y++) {
				SetPixel(hdc, ax, y, col);
			}
		}
		else {
			for (int y = ay; y <= by; y++) {
				int x = (double)(y - ay) * (bx - ax) / (by - ay) + ax;
				SetPixel(hdc, x, y, col);
			}
		}
	}
	else { // Primary axis X
		if (dx < 0) swap(a, b);
		int ax = a.first, bx = b.first;
		int ay = a.second, by = b.second;

		if (dy == 0) {
			for (int x = ax; x <= bx; x++) {
				SetPixel(hdc, x, ay, col);
			}
		}

		else {
			for (int x = ax; x <= bx; x++) {
				int y = (double)(x - ax) * (by - ay) / (bx - ax) + ay;
				SetPixel(hdc, x, y, col);
			}
		}
	}
}

vector<Point> plotGraph(Polinom &f, double hi = 1, double lo = 0, double scale = 100) {

	vector<Point> res;

	for (int i = lo * scale; i < hi * scale - 1; i++) {
		double x  = i / scale;
		double fx = f(x);

		res.push_back(Point(x, fx));

	}

	return res;
}
vector<Point> plotGraph(double (*f)(double), double hi = 1, double lo = 0, double scale = 100) {

	vector<Point> res;

	for (int i = lo * scale; i < hi * scale - 1; i++) {
		double x  = i / scale;
		double fx = f(x);

		res.push_back(Point(x, fx));

	}

	return res;
}


void drawAxis(HDC &hdc, int zerox, int zeroy, int scale) {
	for (int i = zerox - scale * 1.5; i < zerox + scale * 1.5; i++) {

		if (i < 0) continue;

		SetPixel(hdc, i, zeroy, RGB(255, 255, 255));
		if (i == zerox + scale || i == zerox - scale) {
			SetPixel(hdc, i, zeroy + 1, RGB(255, 255, 255));
			SetPixel(hdc, i, zeroy - 1, RGB(255, 255, 255));
		}
	}

	for (int i = zeroy - scale * 1.5; i < zeroy + scale * 1.5; i++) {

		if (i < 0) continue;

		SetPixel(hdc, zerox, i, RGB(255, 255, 255));
		if (i == zeroy + scale || i == zeroy - scale) {
			SetPixel(hdc, zerox + 1, i, RGB(255, 255, 255));
			SetPixel(hdc, zerox - 1, i, RGB(255, 255, 255));
		}
	}
}
void drawGraph(HDC &hdc, vector<Point> dots, double scale, COLORREF col, int zx, int zy) {

	drawAxis(hdc, zx, zy, scale);

	for (int i = 0; i < dots.size() - 1; i++) {
		
		int x1 = dots.at(i).first  * scale + zx;
		int y1 = zy - dots.at(i).second * scale;

		int x2 = dots.at(i + 1).first  * scale + zx;
		int y2 = zy - dots.at(i + 1).second * scale;
	
		drawLine(hdc, Point(x1, y1), Point(x2, y2), col);
	}
}
#undef Point

#endif
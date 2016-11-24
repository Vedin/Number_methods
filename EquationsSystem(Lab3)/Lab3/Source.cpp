#include<iostream>
#include<algorithm>
#include<vector>
#include<ctime>
using namespace std;
const int PS_RAND = 10;
const int OUT_PRES = 5;
const double EPS = 1e-5;
const int INF = 1000000000;
vector<vector<double> > operator *(vector<vector<double> > a, vector<vector<double> > b) {
	vector<vector<double> > ans(a.size(), vector<double> (b[0].size(), 0));
	for (int i = 0; i<a.size(); i++) {
		for (int l = 0; l<b[0].size(); l++) {
			double s = 0;
			for (int j = 0; j<a[i].size(); j++) {
				s += a[i][j] * b[j][l];
			}
			ans[i][l] = s;
		}
	}
	return ans;
}
class EqSys
{
	vector<vector<double> > sys;
	vector<double> B;
	vector<double> x;
	const int RAND_ADD = 10;
	const int MIT = 10000;

public:
	EqSys(int n, int type);
	double Residual(vector<double> rez) {
		vector<vector<double> > prRez(rez.size(), vector<double>(1, 0));
		for (int i = 0; i < prRez.size(); i++) {
			prRez[i][0] = rez[i];
		}
		vector<vector<double> > Ax = sys*prRez;
		vector<double> kAx(Ax.size());
		for (int i = 0; i < Ax.size(); i++) {
			kAx[i] = Ax[i][0];
		}
		return norm(kAx, B);
	}
	
	vector<double> genx0() {
	vector<double> x0(x.size());
	for (int i = 0; i < x.size(); i++) {
		x0[i] = x[i] + (double(rand()) / RAND_MAX) + rand() % RAND_ADD*(rand() % 2 ? 1 : -1);
	}
	return x0;
}

	void out() {
		cout << "A|b:" << endl;
		for (int i = 0; i < sys.size(); i++) {
			for (int j = 0; j < sys[i].size(); j++) {
				cout.precision(OUT_PRES);
				cout.width(OUT_PRES*2+1);
				cout <<fixed<< sys[i][j];
			}
			cout << " | ";
			cout.precision(OUT_PRES);
			cout.width(OUT_PRES*2+1);
			cout << fixed << B[i]<<endl;
		}
		cout << "x:" << endl;
		for (int i = 0; i < x.size(); i++) {
			cout.precision(OUT_PRES);
			cout.width(OUT_PRES * 2 + 1);
			cout << fixed << x[i];
		}
		cout << endl;
	}
	vector<double> tridiagonalEqSolve() {
		vector<double> a(1,0), b, c, d(B);
		for (int i = 0; i < sys.size(); i++) {
			b.push_back(sys[i][i]);
			if (i > 0) a.push_back(sys[i][i - 1]);
			if (i < sys.size() - 1) c.push_back(sys[i][i + 1]);
		}
		c.push_back(0);
		return TridiagonalFull(a, b, c, d, sys.size());
	}
private:
	double norm(const vector<double> &xn,const vector<double> &x0) {
		double ans= 0;
		for (int i = 0; i < xn.size(); i++) {
			ans += sqrt((xn[i] - x0[i])*(xn[i] - x0[i]));
		}
		return ans;
	}
	bool checkEps(vector<double> xn, vector<double> x0) {
		return norm(xn,x0) < EPS;
	}
	double double_rand() {
		return double(rand()) / double(RAND_MAX) + double(rand() % PS_RAND) * (rand() % 2 ? 1 : -1);
	}
	vector<double> TridiagonalFull(const vector<double> a, const vector<double> b, vector<double> c, vector<double> d, int n) {

		vector<double> ans(sys.size());
		c[0] /= b[0];	/* Division by zero risk. */
		d[0] /= b[0];	/* Division by zero would imply a singular matrix. */
		for (unsigned int i = 1; i < n; i++) {
			double id = 1.0 / (b[i] - c[i - 1] * a[i]);  /* Division by zero risk. */
			c[i] *= id;	                         /* Last value calculated is redundant. */
			d[i] = (d[i] - d[i - 1] * a[i]) * id;
		}

		/* Now back substitute. */
		ans[n - 1] = d[n - 1];
		for (int i = n - 2; i >= 0; i--)
			ans[i] = d[i] - c[i] * ans[i + 1];
		return ans;
	}
	void randc(int n) {
		if (n < 2) return;
		sys[0][0] = double_rand();
		sys[0][1] = double_rand();
		sys[n - 1][n - 1] = double_rand();
		sys[n - 1][n - 2] = double_rand();
		for (int i = 1; i < n-1; i++) {
			double sum = 0;
			for (int j = 0; j < 3; j++) {
				sys[i][i+j-1] = double_rand();
				//sum += fabs(sys[i][j]);
			}
			//sys[i][i] = sum;
		}
	}
	/*void gilb(int n) {
		for (int i = 0; i < n; i++) {
			double sum = 0;
			for (int j = 0; j < n; j++) {
				sys[i][j] = 1.0 / double(i + j + 1);
				sum += fabs(sys[i][j]);
			}
			//sys[i][i] = sum;
		}
	}*/
	void genans() {
		for (int i = 0; i < x.size(); i++) {
			x[i] = rand() % PS_RAND;
		}
	}
	void genb() {
		
		vector<vector<double> > prx(x.size(), vector<double> (1,0));
		for (int i = 0; i < prx.size(); i++) {
			prx[i][0] = x[i];
		}
		vector<vector<double> > prb = sys*prx;
		for (int i = 0; i < prb.size(); i++) {
			B[i] = prb[i][0];
		}
	}
};

void outVd(vector<double> a) {
	for (int i = 0; i < a.size(); i++) {
		cout.precision(OUT_PRES);
		cout.width(OUT_PRES * 2 + 1);
		cout << fixed << a[i]<<" ";
	}
	cout << endl;
}

int main() {
	
	cout << "Enter number of equations" << endl;
	int n;
	cin >> n;
	//cout << "Enter type (1 for random 2 for Gilbert matrix)" << endl;
	int type=1;
	//cin >> type;
	freopen("output.txt", "w", stdout);
	srand(time(NULL));
	EqSys A(n, type);
	A.out();
	cout << "Progonka: " << endl;
	vector<double> ans=A.tridiagonalEqSolve();
	outVd(ans);
	cout << "Resudial: " << A.Residual(ans) << endl;
	return 0;
}

EqSys::EqSys(int n, int type)
{
	sys.resize(n, vector<double>(n));
	x.resize(n);
	B.resize(n);
	if (type == 1) randc(n);
	//if (type == 2) gilb(n);
	genans();
	genb();
}

#include<iostream>
#include<vector>
#include<ctime>
using namespace std;
const int PS_RAND = 100;
const double EPS = 1e-4;
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
	vector<double> b;
	vector<double> x;
	const int RAND_ADD = 10;
	const int MIT = 1000;
public:
	EqSys (int n, int type);
	vector<double> Gauss() {
		vector<vector<double> > a=sys;
		for (int i = 0; i < a.size(); i++)
			a[i].push_back(b[i]);
		vector<double> ans(a.size());
		int st=gaussAll(a, ans);
		if (st == INF) {
			cout << "Unlimited system";
		}
		if (st == 0) {
			cout << "No solutions";
		}
		return ans;
	}
	vector<double> genx0() {
		vector<double> x0(x.size());
		for (int i = 0; i < x.size(); i++) {
			x0[i] = x[i] + (double(rand()) / RAND_MAX) + rand() % RAND_ADD*(rand() % 2 ? 1 : -1);
		}
		return x0;
	}
	vector<double> Zeydel() {
		vector<double> x0=genx0();
		vector<double> xn = x0;
		int it = 0;
		do {
			for (int i = 0; i < sys.size(); i++) {
				double sum = 0;
				for (int j = 0; j < sys[i].size(); j++) {
					if (i != j) {
						sum = sum + x0[j] * sys[i][j];
					}
				}
				x0[i] = (b[i] - sum) / sys[i][i];
			}
			swap(x0, xn);
			it++;
		} while (!checkEps(x0, xn) && it<MIT);
		return x0;
	}
	vector<double> Yakobi() {
		vector<double> x0=genx0();
		vector<double> xn = x0;
		int it = 0;
		do {
			for (int i = 0; i < sys.size(); i++) {
				double sum = 0;
				for (int j = 0; j < sys[i].size(); j++) {
					if (i != j) {
						sum = sum + x0[j] * sys[i][j];
					}
				}
				xn[i] = (b[i] - sum) / sys[i][i];
			}
			swap(x0, xn);
			it++;
		} while (!checkEps(x0, xn) && it<MIT);
		return x0;
	}
private:
	bool checkEps(vector<double> xn, vector<double> x0) {
		double norm = 0;
		for (int i = 0; i < xn.size(); i++) {
			norm += sqrt((xn[i] - x0[i])*(xn[i] - x0[i]));
		}
		return norm < EPS;
	}
	int gaussAll(vector < vector<double> > a, vector<double> & ans) {
		int n = (int)a.size();
		int m = (int)a[0].size() - 1;

		vector<int> where(m, -1);
		for (int col = 0, row = 0; col<m && row<n; ++col) {
			int sel = row;
			for (int i = row; i<n; ++i)
				if (abs(a[i][col]) > abs(a[sel][col]))
					sel = i;
			if (abs(a[sel][col]) < EPS)
				continue;
			for (int i = col; i <= m; ++i)
				swap(a[sel][i], a[row][i]);
			where[col] = row;

			for (int i = 0; i<n; ++i)
				if (i != row) {
					double c = a[i][col] / a[row][col];
					for (int j = col; j <= m; ++j)
						a[i][j] -= a[row][j] * c;
				}
			++row;
		}

		ans.assign(m, 0);
		for (int i = 0; i<m; ++i)
			if (where[i] != -1)
				ans[i] = a[where[i]][m] / a[where[i]][i];
		for (int i = 0; i<n; ++i) {
			double sum = 0;
			for (int j = 0; j<m; ++j)
				sum += ans[j] * a[i][j];
			if (abs(sum - a[i][m]) > EPS)
				return 0;
		}

		for (int i = 0; i<m; ++i)
			if (where[i] == -1)
				return INF;
		return 1;
	}
	void randc(int n) {
		for (int i = 0; i < n; i++) {
			double sum = 0;
			for (int j = 0; j < n; j++) {
				sys[i][j] = double(rand()) / double(RAND_MAX)+double(rand()%PS_RAND) * (rand()%2 ? 1 : -1);
				sum += sys[i][j];
			}
			sys[i][i] = sum;
		}
	}
	void gilb(int n) {
		for (int i = 0; i < n; i++) {
			double sum = 0;
			for (int j = 0; j < n; j++) {
				sys[i][j] = 1.0 / double(i + j + 2);
				sum += sys[i][j];
			}
			sys[i][i] = sum;
		}
	}
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
			b[i] = prb[i][0];
		}
	}
};


int main() {
	srand(time(NULL));
	EqSys A(2, 1);
	vector<double> ansG = A.Gauss();
	vector<double> ansYa = A.Yakobi();
	vector<double> ansZe = A.Zeydel();
	return 0;
}

EqSys::EqSys(int n, int type)
{
	sys.resize(n, vector<double>(n));
	x.resize(n);
	b.resize(n);
	if (type == 1) randc(n);
	if (type == 2) gilb(n);
	genans();
	genb();
}

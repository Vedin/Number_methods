#include<iostream>
#include<vector>
#include<cmath>
#include<iomanip>
#include<cstdio>
#define INF 1000010000
#define MAXN 5050
#define alpha 0.85
#define eps 1e-7
using namespace std;

int n, m; 
vector<int> g[MAXN];
vector<vector<double> > H(MAXN, vector<double>(MAXN, 0.));
vector<bool> a(MAXN, 0);
vector<vector<double> > G(MAXN, vector<double>(MAXN, 0.));

inline double EuclideanNorm(vector<double> v) {
	double sum = 0;
	for (int i = 0; i < v.size(); i++) sum += v[i] * v[i];
	return sqrt(sum);
}

inline double l_Norm(vector<double> v) {
	double sum = 0;
	for(int i=0; i<n; i++)sum += abs(v[i]);
	return sum;
}

inline vector<double> Substract(vector<double> & a, vector<double> & b) {
	vector<double> res(n, 0.);
	for(int i=0; i<n; i++)res[i] = a[i] - b[i];
	return res;
}

inline void Build_H() {
	for(int i=0; i<n; i++) {
		a[i] = !(g[i].size());
		double cur = 1. / (double)(g[i].size());
		for (int j : g[i]) {
			H[i][j] = cur;
		}
	}
}

inline void Build_G() {
	for(int i=0; i<n; i++) {
		double t = (alpha * a[i] + (1 - alpha)) / (double)n;
		for(int j=0; j<n; j++) {
			G[i][j] = alpha * H[i][j];
			G[i][j] += t;
		}
	}
}

inline void Transpose(vector<vector<double> > & A) {
	for(int i=0; i<n; i++) {
		for (int j = i + 1; j < n; j++) {
			swap(A[i][j], A[j][i]);
		}
	}
}

inline vector<double> Gauss(vector<vector<double> > A, vector<double> B) {//Gauss
	vector<double> X(n);
	vector<bool> isZero(n, 0);
	for(int j=0; j<n; j++) {
		if (isZero[j])continue;
		double k = A[j][j];
		for(int i=0; i<n; i++) {
			A[j][i] /= k;
		}
		B[j] /= k;
		for (int i = j + 1; i < n; i++) {
			double k = A[i][j] / A[j][j];
			for (int r = 0; r < n; r++) {
				A[i][r] -= k * A[j][r];
			}

			B[i] -= k * B[j];
		}
		for (int ii = 0; ii<n; ii++) {
			for (int jj = 0; jj<n; jj++) {
				if (abs(A[ii][jj]) > eps) { isZero[ii] = 0; break; }
				else isZero[ii] = 1;
			}
		}
	}
	for (int i = n - 1; i >= 0; i--) {
		if (isZero[i]) {
			X[i] = 1.;
			continue;
		}
		double res = B[i];
		for (int j = i + 1; j < n;j++) {
			res -= A[i][j] * X[j];
		}
		res /= A[i][i];
		X[i] = res;
		if (X[i] == -0.)X[i] = 0;
	}

	return X;
}

inline vector<double> MultV(vector<vector<double> > A, vector<double> v) {
	vector<double> ans((v.size()), 0);
	for (int i = 0; i < (v.size()); i++) {
		for (int j = 0; j < v.size(); j++) {
			ans[i] += A[i][j] * v[j];
		}
	}
	return ans;
}

inline vector<double> Step(vector<vector<double> > G) {
	vector<double> X(n, 1);
	vector<double> Y(n, 1);
	int iter = 0;
	double cur_lambda1, cur_lambda2 = INF;
	do {
		iter++;
		X = Y;
		cur_lambda1 = cur_lambda2;
		Y = MultV(G, Y);
		int pos = 0;
		while (abs(X[pos]) < eps)++pos;
		cur_lambda2 = Y[pos] / X[pos];
		double norm = l_Norm(Y);
		for(int i=0; i<n; i++)Y[i] /= norm;

	} while (abs(cur_lambda1 - cur_lambda2) > eps);
	cout << "Iter: " << iter << endl;
	return X;
}

inline vector<double> PageRank(vector<vector<double> > G) {
	Transpose(G);
	vector<double> v = Step(G);
	return v;
}

int main() {
	freopen("input.txt", "r", stdin);
	cout.precision(3);
	cout.setf(ios::fixed);
	cin >> n >> m;
	for (int i = 0; i < m; i++) {
		int u, v;
		cin >> u >> v;
		--u;
		--v;
		g[u].push_back(v);
	}
	Build_H();
	Build_G();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << setw(4) << G[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	vector<double> v = PageRank(G);
	for (int i = 0; i < n;i++) cout << setw(4) << v[i] << " ";

	return 0;
}

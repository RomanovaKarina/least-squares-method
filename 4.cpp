#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

const int m = 3;
double* methodGauss(double[m + 1][m + 1], double*);

int main() {
	const int N=21;

	double t[] = {
		0, 5,10,15,20,
		25,30,35,40,45,
		50,55,60,65,70,
		75,80,85,90,95,100
	};
	double C[] = {
		1.00762, 1.00392, 1.00153, 1 , 0.99907, 0.99852,
		0.99826, 0.99818, 0.99828, 0.99849, 0.99878,
		0.999919, 0.99967, 1.00024, 1.00091, 1.00167, 1.00253,
		1.00351,1.00461, 1.00586, 1.00721
	};

	double powerX[2 * m] = { 0 }, sumX[m + 1][m + 1] = { 0 }, praw[m + 1] = { 0 }, S2 = 0;

	cout <<"t:" << endl;
	for (int i = 0; i < N; i++) {
		cout  << t[i] << setw(5);
		if (i % 5 == 0 && i!=0)
			cout << endl;
	}
	cout << endl;
	cout << "C:" << endl;
	for (int i = 0; i < N; i++) {
		cout << C[i] << "  ";
		if (i % 5 == 0 && i != 0)
			cout << endl;
	}
	cout << endl;
	cout << "POWERX:" << endl;
	for (int k = 0; k < 2 * m; k++)
	{
		powerX[k] = 0;
		for (int i = 0; i < N; i++)
		{
			powerX[k] += pow(t[i], k + 1);
		}
		cout << powerX[k] << " ";
	}
	cout << endl;
	cout << endl;

	for (int i = 0; i < m + 1; i++)
		for (int j = 0; j < m + 1; j++)
			sumX[i][j] = powerX[i + j ];

	sumX[0][0] = N;

	cout << "SUMX:" << endl;
	for (int i = 0; i < m + 1; i++) {
		for (int j = 0; j < m + 1; j++)
			cout << sumX[i][j] << " ";
		cout << endl;
	}
	cout << endl;

	for (int i = 0; i < m + 1; i++)
	{
		praw[i] = 0;
		for (int j = 0; j < N; j++)
		{
			praw[i] += C[j] * pow(t[j], i);
		}
	}

	cout << "PRAW:" << endl;
	for (int i = 0; i < m+1; i++)
	{
		cout << praw[i] << " ";
	}
	cout << endl;

	double* a = methodGauss(sumX, praw);

	cout << "a:" << endl;
	for (int i = 0; i < m + 1; i++) {
		cout << a[i] << " ";
	}
	cout << endl;

	for (int i = 0; i < N; i++)
	{
		double sum = C[i];
		for (int j = 0; j < m + 1; j++)
		{
			sum -= a[j] * pow(t[i], j);
		}
		S2 += pow(sum, 2);
	}
	S2 /= N - m - 1;
	double sigma = sqrt(S2);
	return 0;
}

double* methodGauss(double mA[m + 1][m + 1], double* cB)
{
	const int n = m + 1;
	double X[n] = { 0 };
	for (int i = 0; i < n; i++) {
		int maxIndex = i;
		double max = mA[i][i];
		for (int j = i + 1; j < n; j++) {
			if (abs(max) < abs(mA[j][i])) {
				maxIndex = j;
				max = mA[j][i];
			}
		}
		if (i != maxIndex) {
			double root = cB[i];
			cB[i] = cB[maxIndex];
			cB[maxIndex] = root;
			for (int j = 0; j < n; j++) {
				double x = mA[i][j];
				mA[i][j] = mA[maxIndex][j];
				mA[maxIndex][j] = x;
			}
		}
		double a = mA[i][i];
		for (int j = i; j < n; j++)
		{
			mA[i][j] /= a;
		}
		cB[i] /= a;
		for (int j = i + 1; j < n; j++)
		{
			double s = mA[j][i];
			for (int k = i; k < n; k++)
			{
				mA[j][k] -= s * mA[i][k];
			}
			cB[j] -= s * cB[i];
		}
	}

	for (int k = n - 1; k >= 0; k--)
	{
		X[k] = cB[k];
		for (int i = n - 1; i > k; i--)
		{
			X[k] -= mA[k][i] * X[i];
		}
	}
	cout << endl;
	return X;
}
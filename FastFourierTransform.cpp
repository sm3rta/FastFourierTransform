#include <iostream>
#include <iomanip>
#include <complex>
#include <complex.h>
using namespace std;

#define PI 3.14159265358979323846

void print(complex<double>n) {
	if (n.imag() < 0)
		cout << n.real() << "-j" << 0 - n.imag();
	else
		cout << n.real() << "+j" << n.imag();
}

complex<double> w_n_k(int N, int k) {
	return complex<double>(cos(2 * PI*k / N), -sin(2 * PI*k / N));
}

complex<double>* FFT(complex <double>* arr, int N) {
	if (N == 2) {
		complex <double> *x = new complex <double>[2];

		x[0] = arr[0] + w_n_k(2, 0)*arr[1];
		x[1] = arr[0] - w_n_k(2, 0)*arr[1];

		return x;
	}

	complex <double>* odd = new complex <double>[N / 2];
	complex <double>* even = new complex <double>[N / 2];

	for (int i = 0; i < N / 2; i++) {
		even[i] = arr[(2 * i)];
		odd[i] = arr[(2 * i + 1)];
	}

	complex <double>* Odd = FFT(odd, N / 2);
	complex <double>* Even = FFT(even, N / 2);


	complex <double>* temp = new complex<double>[N];

	for (int i = 0; i < N; i++) {
		if (i < N / 2)
			temp[i] = Even[i] + Odd[i] * w_n_k(N, i);
		else {
			int k = i - (N / 2);
			temp[i] = Even[k] - Odd[k] * w_n_k(N, k);
		}
	}
	return temp;
}

void main() {
	cout << fixed << setprecision(1);

	int N;

	cout << "Enter N: ";

	while (cin >> N) {
		if (log2(N) == int(log2(N)))
			break;
		else
			cout << "ERROR: log2(N) must be an integer" << endl;
	}

	complex<double> * x = new complex<double>[N];
	cout << "Enter x(n): ";
	int j;
	for (int i = 0; i < N; i++) {
		cin >> j;
		x[i].real(j);
		x[i].imag(0);
	}

	complex<double> *X = FFT(x, N);

	//print X(k)
	cout << "\nX(k) = \t\t";
	for (int i = 0; i < N; i++) {
		print(X[i]);
		cout << "\t";
	}
	cout << endl;

	//print |X(k)|
	cout << "|X(k)| = \t";
	for (int i = 0; i < N; i++)
		cout << cabs(X[i]) << "\t\t";
	cout << endl;

	//print arg(X(k))
	cout << "arg(X(k)) = \t";
	for (int i = 0; i < N; i++)
		cout << carg(X[i]) << "\t\t";
	cout << endl;
}

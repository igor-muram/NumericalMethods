#include <cmath>
#include <cstdio>

double rho = 1000.0;
double Csnd = 1260.0;
double l = 1.0;
double d = 0.01;
double alpha = 3.14159265359 / 4;
double q0 = 1e-3;
double S = alpha * d * d;
double V = l * S;
double Patm = 1e+5;
double Es = rho * Csnd * Csnd;
double C = V / Es;
double zeta = 1 - cos(alpha);
double B = 1 / l * sqrt(zeta / (2 * rho));
double F = S * sqrt(2 / rho);

double Qpump(double t)
{
	return t < 1.0 ? q0 * t : q0;
}

double dp(double dq)
{
	return dq / C;
}

double dq(double q, double dp)
{
	double signDP = dp > 0 ? 1 : -1;

	return B * sqrt(abs(dp)) * (F * sqrt(abs(dp) / zeta) * signDP - q);
}

void Runge1(double a, double b, double h, double p0, double q0)
{
	int n = (b - a) / h;
	double K[4][2] = { 0 };

	for (int i = 0; i < n; i++)
	{
		double t = i * h;

		K[0][0] = dp(Qpump(t) - q0);
		K[1][0] = dp(Qpump(t + h / 2) - (q0));
		K[2][0] = dp(Qpump(t + h / 2) - (q0));
		K[3][0] = dp(Qpump(t + h) - (q0));
		p0 += h * (K[0][0] + 2 * K[1][0] + 2 * K[2][0] + K[3][0]) / 6;

		K[0][1] = dq(q0, p0 - Patm);
		K[1][1] = dq(q0 + h * K[0][1] / 2, p0 - Patm);
		K[2][1] = dq(q0 + h * K[1][1] / 2, p0 - Patm);
		K[3][1] = dq(q0 + h * K[2][1], p0 - Patm);
		q0 += h * (K[0][1] + 2 * K[1][1] + 2 * K[2][1] + K[3][1]) / 6;
	}

	printf_s("p: %.5f\t\tq: %.5f\n", p0, q0);
}

void Runge2(double a, double b, double h, double p0, double q0)
{
	int n = (b - a) / h;
	double K[4][2] = { 0 };

	for (int i = 0; i < n; i++)
	{
		double t = i * h;

		K[0][0] = dp(Qpump(t) - q0);
		K[0][1] = dq(q0, p0 - Patm);

		K[1][0] = dp(Qpump(t + h / 2) - (q0 + h * K[0][1] / 2));
		K[1][1] = dq(q0 + h * K[0][1] / 2, p0 + h * K[0][0] / 2 - Patm);

		K[2][0] = dp(Qpump(t + h / 2) - (q0 + h * K[1][1] / 2));
		K[2][1] = dq(q0 + h * K[1][1] / 2, p0 + h * K[1][0] / 2 - Patm);

		K[3][0] = dp(Qpump(t + h) - (q0 + h * K[2][1]));
		K[3][1] = dq(q0 + h * K[2][1], p0 + h * K[2][0] - Patm);

		p0 += h * (K[0][0] + 2 * K[1][0] + 2 * K[2][0] + K[3][0]) / 6;
		q0 += h * (K[0][1] + 2 * K[1][1] + 2 * K[2][1] + K[3][1]) / 6;
	}

	printf_s("p: %.5f\t\tq: %.5f\n", p0, q0);
}

int main()
{
	Runge1(0.0, 20.0, 1e-5, Patm, 0.0);
	Runge2(0.0, 20.0, 1e-5, Patm, 0.0);
	return 0;
}
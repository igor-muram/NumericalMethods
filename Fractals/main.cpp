#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <string>
#include <complex>

#include "GL\freeglut.h"

using namespace std;

int n = 1, key = 0;
double x2, y2, angle = 0.5;	// for tree
int r = 2;					// for tree and Pifagor

int iteration = 0;

int pifagorR = 1;
double q = 1.0;

vector<double> x, y;

enum {
	KochCurveKey,
	KochSnowflakeKey,
	KochAntiSnowflakeKey,
	TreeKey,
	LevyCurveKey,
	LevySquare1Key,
	LevySquare2Key,
	PifagorKey,
	DragonLevyCurveKey,
	DragonLevySquare1Key,
	DragonLevySquare2Key,
	MandelbrotKey
};

void KochPattern(double x1, double y1, double x2, double y2, int n, bool f)
{
	x[n] = x1;
	y[n] = y1;

	x[n + 1] = (2 * x1 + x2) / 3.0;
	y[n + 1] = (2 * y1 + y2) / 3.0;

	x[n + 2] = (x2 + x1) / 2.0 + (f ? 1 : -1) * (y1 - y2) / (2 * sqrt(3.0));
	y[n + 2] = (y2 + y1) / 2.0 + (f ? 1 : -1) * (x2 - x1) / (2 * sqrt(3.0));

	x[n + 3] = (x1 + 2 * x2) / 3.0;
	y[n + 3] = (y1 + 2 * y2) / 3.0;

	x[n + 4] = x2;
	y[n + 4] = y2;
}

void KochCalculate(int size, int k, int coeff, bool f)
{
	for (int j = 0; j <= k; j++)
	{
		vector<double> tempX, tempY;
		tempX.resize(size);
		tempY.resize(size);

		n = coeff * pow(4.0, j + 1);

		for (int i = 0; i <= n; i++)
		{
			tempX[i] = x[i];
			tempY[i] = y[i];
		}

		for (int i = 0; i < n; i++)
			KochPattern(tempX[i], tempY[i], tempX[i + 1], tempY[i + 1], 4 * i, f);
	}
}

void KochCurve(int count, bool f)
{
	int size = pow(4.0, count + 2) + 1;
	x.resize(size);
	y.resize(size);

	KochPattern(-1.0, 0.0, 1.0, 0.0, 0, f);

	KochCalculate(size, count, 1, f);
}

void KochSnowflake(int count, bool f)
{
	int size = 3 * pow(4.0, count + 2) + 1;
	x.resize(size);
	y.resize(size);

	double Ax = -0.7, Ay = -0.4, Bx = 0.7, By = -0.4, Cx, Cy;

	Cx = (Bx + Ax) / 2 + (Ay - By) * sqrt(3.0) / 2;
	Cy = (By + Ay) / 2 + (Bx - Ax) * sqrt(3.0) / 2;

	KochPattern(Bx, By, Ax, Ay, 0, f);
	KochPattern(Ax, Ay, Cx, Cy, 4, f);
	KochPattern(Cx, Cy, Bx, By, 8, f);

	KochCalculate(size, count, 3, f);
}

void TreeCalculate(double x1, double y1, double a, double L, int count)
{
	for (int i = 0; i <= 2 * r; i += r)
	{
		x2 = x1 + cos(a + i * angle) * L;
		y2 = y1 + sin(a + i * angle) * L;

		glVertex2d(x1, y1);
		glVertex2d(x2, y2);

		if (count > 0)
			TreeCalculate(x2, y2, a + (i - 1) * angle, L / 2, count - 1);
	}
}

void LevyCalculate(double x1, double y1, double x2, double y2, int count, bool f)
{
	if (count == 0)
	{
		glVertex2d(x1, y1);
		glVertex2d(x2, y2);
	}
	else
	{
		double x3 = (x1 + x2) / 2.0 + (y2 - y1) / 2.0;
		double y3 = (y1 + y2) / 2.0 - (x2 - x1) / 2.0;
		LevyCalculate(x1, y1, x3, y3, count - 1, f);

		if (f)
			LevyCalculate(x3, y3, x2, y2, count - 1, f);
		else
			LevyCalculate(x2, y2, x3, y3, count - 1, f);
	}
}

void PifagorPattern(double x1, double y1, double x2, double y2, int count)
{
	x[count + 4] = x1;
	y[count + 4] = y1;
	x[count + 3] = x2;
	y[count + 3] = y2;

	double x3 = x[count] = x1 + y1 - y2;
	double y3 = y[count] = y1 + x2 - x1;
	double x4 = x[count + 2] = x2 + y1 - y2;
	double y4 = y[count + 2] = y2 + x2 - x1;

	x[count + 1] = (x3 + x4 * q * q + (y3 - y4) * q) / (1 + q * q);
	y[count + 1] = (y3 + y4 * q * q + (x4 - x3) * q) / (1 + q * q);
}

void PifagorCalculate(int j)
{
	int count = 5 * (pow(2, j) - 1);
	int size = 10 * (count + 1);
	x.resize(size);
	y.resize(size);

	PifagorPattern(-0.1, -0.5, 0.15, -0.5, 0);

	for (int i = 0; i < count; i += 5)
	{
		PifagorPattern(x[i], y[i], x[i + 1], y[i + 1], 5 * pifagorR);
		PifagorPattern(x[i + 1], y[i + 1], x[i + 2], y[i + 2], 5 * (pifagorR + 1));
		pifagorR += 2;
	}
}

void Reshape(int w, int h)
{
	glViewport(0, 0, w, h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, w, 0, h);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void MandelbrotDraw()
{
	glClearColor(0.11, 0.22, 0.31, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);

	glBegin(GL_POINTS);
	for (int y = 0; y < 800; y++)
		for (int x = 0; x < 1200; x++)
		{
			complex<double> z(0, 0);
			int i = 0;
			for (; i < 500 && abs(z) < 2; i++)
				z = z * z + complex<double>((x - 750) / 300.0, (y - 400) / 300.0);

			double r = 0.1 + i * 0.03 * 0.2;
			double g = 0.2 + i * 0.03 * 0.3;
			double b = 0.3 + i * 0.03 * 0.1;

			glColor3d(r, g, b);
			glVertex2d(x, y);
		}
	glEnd();

	glutSwapBuffers();
}

void Draw()
{
	glClearColor(0.11, 0.22, 0.31, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3d(1, 1, 1);
	glLineWidth(2.0);

	glBegin(GL_LINES);
	switch (key - 1)
	{
	case KochCurveKey:
	case KochSnowflakeKey:
	case KochAntiSnowflakeKey:
	{
		switch (key - 1)
		{
		case KochCurveKey:
			KochCurve(iteration, true);
			break;

		case KochSnowflakeKey:
			KochSnowflake(iteration, true);
			break;

		case KochAntiSnowflakeKey:
			KochSnowflake(iteration, false);
			break;

		}

		for (int i = 0; i < 4 * n; i++)
		{
			glVertex2d(x[i], y[i]);
			glVertex2d(x[i + 1], y[i + 1]);
		}
	}
	break;

	case TreeKey:
	{
		glVertex2d(0.0, -1.0);
		glVertex2d(0.0, -0.1);
		TreeCalculate(0.0, -0.1, 0.60, 0.5, iteration);
	}
	break;

	case LevyCurveKey:
	{
		LevyCalculate(0.4, 0.0, -0.4, 0.0, iteration, true);
	}
	break;

	case LevySquare1Key:
	case DragonLevySquare1Key:
	{
		bool flag = (key == 6);

		LevyCalculate(0.4, 0.0, -0.4, 0.0, iteration, flag);
		LevyCalculate(-0.4, 0.0, 0.4, 0.0, iteration, flag);
	}
	break;

	case LevySquare2Key:
	case DragonLevySquare2Key:
	{
		bool flag = (key == 7);

		LevyCalculate(0.3, 0.3, -0.3, 0.3, iteration, flag);
		LevyCalculate(-0.3, -0.3, 0.3, -0.3, iteration, flag);
		LevyCalculate(-0.3, 0.3, -0.3, -0.3, iteration, flag);
		LevyCalculate(0.3, -0.3, 0.3, 0.3, iteration, flag);
	}
	break;

	case DragonLevyCurveKey:
	{
		LevyCalculate(0.5, 0.0, -0.7, -0.3, iteration, false);
	}
	break;

	case PifagorKey:
	{
		PifagorCalculate(iteration);

		glVertex2d(x[3], y[3]);
		glVertex2d(x[4], y[4]);

		for (int i = 0; i < 5 * pifagorR; i += 5)
		{
			glVertex2d(x[i], y[i]);
			glVertex2d(x[i + 1], y[i + 1]);

			glVertex2d(x[i], y[i]);
			glVertex2d(x[i + 2], y[i + 2]);

			glVertex2d(x[i + 1], y[i + 1]);
			glVertex2d(x[i + 2], y[i + 2]);

			glVertex2d(x[i + 2], y[i + 2]);
			glVertex2d(x[i + 3], y[i + 3]);

			glVertex2d(x[i + 4], y[i + 4]);
			glVertex2d(x[i], y[i]);
		}
	}
	break;
	}
	glEnd();
	glFinish();
}

void Keyboard(unsigned char k, int x, int y)
{
	if (k == 'w')
		iteration++;

	if (k == 's')
	{
		if (iteration > 0)
			iteration--;

		pifagorR = 1;
	}

	if (k == 'a')
		if (key > 1)
		{
			key--;
			iteration = 0;
		}

	if (k == 'd')
		if (key < 11)
		{
			key++;
			iteration = 0;
		}

	if (k == 'r')
	{
		if (r < 12)
			r++;
	}

	if (k == 'f')
		if (r > 1)
			r--;

	if (k == 't')
		angle += 0.05;

	if (k == 'g')
		angle -= 0.05;

	if (k == 'y')
	{
		q += 0.01;
		pifagorR = 1;
	}

	if (k == 'h')
	{
		if (q > 0.1)
			q -= 0.01;
		else
			q = 0.1;

		pifagorR = 1;
	}

	glutPostRedisplay();
}

void InitWindow(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);
	if (key == 12)
		glutInitWindowSize(1200, 800);
	else
		glutInitWindowSize(900, 900);
	glutCreateWindow("Fractals");
	glutKeyboardFunc(Keyboard);
}

void main(int argc, char* argv[])
{
	cout << "1. Koch curve" << endl;
	cout << "2. Koch snowflake" << endl;
	cout << "3. Koch antisnowflake" << endl;
	cout << "4. Tree" << endl;
	cout << "5. Levy curve" << endl;
	cout << "6. Levy square 1" << endl;
	cout << "7. Levy square 2" << endl;
	cout << "8. Pifagor tree" << endl;
	cout << "9. Dragon Levy curve" << endl;
	cout << "10. Dragon Levy square 1" << endl;
	cout << "11. Dragon Levy square 2" << endl;
	cout << "12. Mandelbrot set" << endl;
	cout << "Choose a fractal: ";

	cin >> key;

	while (key < 1 || key > 12)
	{
		cout << endl << "Incorrect value. Try again. Choose a fractal: ";
		cin >> key;
	}

	if (key == 12)
	{
		InitWindow(argc, argv);
		glutReshapeFunc(Reshape);
		glutDisplayFunc(MandelbrotDraw);
		glutMainLoop();
	}
	else
	{
		InitWindow(argc, argv);
		glutDisplayFunc(Draw);
		glutMainLoop();
	}

}
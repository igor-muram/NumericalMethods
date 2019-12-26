#include "GL\freeglut.h"
#include <fstream>
#include <vector>
#include <functional>
#include <string>
using namespace std;

int Width = 1280, Height = 720;
int xAxis = Width / 2, yAxis = Height / 2;
int step = 70;

struct Point
{
	Point(double x, double y) : x(x), y(y) {};
	double x, y;
};

vector<Point> points;

vector<function<double(double, double)>> F =
{
	[](double x, double y) { return (x + 2) * (x + 2) + (y - 1) * (y - 1) - 4; },
	[](double x, double y) { return (x - 1) * (x - 1) + (y - 1) * (y - 1) - 4; },
	[](double x, double y) { return y - x + 3.5; }
};

void ReadPoints(string filename)
{
	ifstream in(filename);
	int n;
	in >> n;
	for (int i = 0; i < n; i++)
	{
		double x, y;
		in >> x >> y;
		points.emplace_back(x, y);
	}

	in.close();
}

double Norm(double x, double y)
{
	double sum = 0;
	for (auto f : F)
		sum += f(x, y) * f(x, y);
	return sqrt(sum);
}

void DrawDifference(double h)
{
	glPointSize(4);
	glBegin(GL_POINTS);
	for (int i = -500; i < 500; i++)
	{
		for (int j = -500; j < 500; j++)
		{
			double x = i * h;
			double y = j * h;

			double diff = Norm(x, y);
			if (diff < 0.5)
			{
				glColor3ub(50, 50, 50);
				glVertex2f(xAxis + x * step, yAxis + y * step);
			}

			if (diff > 0.5 && diff < 2.0)
			{
				glColor3ub(80, 80, 80);
				glVertex2f(xAxis + x * step, yAxis + y * step);
			}

			for (int i = 1; i < 50; i++)
				if (diff > i * 2.0 && diff < (i + 1) * 2.0)
				{
					double color = 80 + i * 5;
					if (color > 255) color = 255;
					glColor3ub(color, color, color);
					glVertex2f(xAxis + x * step, yAxis + y * step);
				}
		}
	}
	glEnd();
}

void DrawAxis()
{
	glColor3ub(0, 0, 0);
	glLineWidth(2.0);
	glBegin(GL_LINES);

	glVertex2f(-10000, yAxis);
	glVertex2f(10000, yAxis);

	for (int i = step; i < 5000; i += step)
	{
		glVertex2f(xAxis - i, yAxis);
		glVertex2f(xAxis - i, yAxis + 3);

		glVertex2f(xAxis + i, yAxis);
		glVertex2f(xAxis + i, yAxis + 3);
	}

	glVertex2f(xAxis, -10000);
	glVertex2f(xAxis, 10000);

	for (int i = step; i < 5000; i += step)
	{
		glVertex2f(xAxis, yAxis - i);
		glVertex2f(xAxis + 3, yAxis - i);

		glVertex2f(xAxis, yAxis + i);
		glVertex2f(xAxis + 3, yAxis + i);
	}
	glEnd();
}

void DrawFunc(double a, double b, double h, function<double(double)> f)
{
	int n = (b - a) / h;

	glLineWidth(2);
	glColor3ub(0, 200, 0);
	glBegin(GL_LINE_STRIP);
	for (int i = 0; i < n; i++)
	{
		double t = a + i * h;
		glVertex2f(xAxis + t * step, yAxis + f(t) * step);
	}
	glEnd();
	glLineWidth(1);
}

void DrawParametric(double a, double b, double h, function<double(double)> xt, function<double(double)> yt)
{
	int n = (b - a) / h;

	glLineWidth(2);
	glColor3ub(0, 200, 0);
	glBegin(GL_LINE_STRIP);
	for (int i = 0; i < n; i++)
	{
		double t = a + i * h;
		glVertex2f(xAxis + xt(t) * step, yAxis + yt(t) * step);
	}
	glEnd();
	glLineWidth(1);
}

void DrawPoints(vector<Point> points)
{
	glLineWidth(2);
	glColor3ub(65, 105, 225);
	glBegin(GL_LINE_STRIP);
	for (auto point : points)
		glVertex2f(xAxis + point.x * step, yAxis + point.y * step);
	glEnd();
	glLineWidth(1);

	glColor3ub(255, 0, 0);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (auto point : points)
		glVertex2f(xAxis + point.x * step, yAxis + point.y * step);
	glEnd();
	glPointSize(1);
}

void Display()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);

	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);

	DrawDifference(0.02);
	DrawAxis();
	DrawParametric(0, 10, 0.1, [](double t) { return -2 + 2 * cos(t); }, [](double t) { return 1 + 2 * sin(t); });
	DrawParametric(0, 10, 0.1, [](double t) { return 1 + 2 * cos(t); }, [](double t) { return 1 + 2 * sin(t); });
	DrawFunc(-10, 10, 0.1, [](double t) { return t - 3.5; });
	//DrawPoints(points);
	glFinish();
}

void Reshape(GLint w, GLint h)
{
	Width = w; Height = h;
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, w, 0, h);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void Keyboard(unsigned char key, int x, int y)
{
	if (key == 'w')
		yAxis -= 25;

	if (key == 's')
		yAxis += 25;

	if (key == 'a')
		xAxis += 25;

	if (key == 'd')
		xAxis -= 25;

	if (key == 'p' && step < 200)
		step += 10;

	if (key == 'm' && step > 10)
		step -= 10;

	if (key == 'x')
		glutExit();

	glutPostRedisplay();
}

void main(int argc, char* argv[])
{
	ReadPoints("C:/input/points.txt");
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);
	glutInitWindowSize(Width, Height);
	glutCreateWindow("Визуализация");
	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutKeyboardFunc(Keyboard);
	glutMainLoop();
}
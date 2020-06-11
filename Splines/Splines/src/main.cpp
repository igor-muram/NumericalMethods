#include "GL\freeglut.h"
#include <vector>
#include <iostream>
#include <string>

#include "Spline.h"
using namespace std;

int Width = 1280, Height = 720;
int xAxis = Width / 2, yAxis = Height / 2;
double step = 200;
int gridWidth = 200;

int currentPoint = -1;
bool drawLabel = false;

vector<Point> points, Bpoints;
vector<Point> spline1, spline2, spline3, spline4, spline5;

bool s1 = false;
bool s2 = false;
bool s3 = false;
bool s4 = false;
bool s5 = false;
bool s6 = true;

void DrawText(string text, int x, int y)
{
	glColor3ub(230, 230, 230);
	glRasterPos2i(x, y);
	for (char c : text)
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, (int)c);
}

void DrawGrid()
{
	glColor3ub(230, 230, 230);
	glLineWidth(2.0);
	glBegin(GL_LINES);

	glVertex2f(-10000, yAxis);
	glVertex2f(10000, yAxis);

	glEnd();

	glLineWidth(1.0);
	glBegin(GL_LINES);
	for (int i = gridWidth; i < 5000; i += gridWidth)
	{
		glVertex2f(xAxis - i, -10000);
		glVertex2f(xAxis - i, 10000);

		glVertex2f(xAxis + i, -10000);
		glVertex2f(xAxis + i, 10000);
	}
	glEnd();


	glLineWidth(2.0);
	glBegin(GL_LINES);

	glVertex2f(xAxis, -10000);
	glVertex2f(xAxis, 10000);

	glEnd();


	glLineWidth(1.0);
	glBegin(GL_LINES);
	for (int i = gridWidth; i < 5000; i += gridWidth)
	{
		glVertex2f(10000, yAxis - i);
		glVertex2f(-10000, yAxis - i);

		glVertex2f(10000, yAxis + i);
		glVertex2f(-10000, yAxis + i);
	}
	glEnd();
}

void DrawLabels()
{
	DrawText("y", xAxis + 10, Height - 10);
	DrawText("x", Width - 10, yAxis + 10);

	for (int i = 0; i < 5000; i += gridWidth)
	{
		DrawText(to_string(-i / step), xAxis - i + 2, 2);
		DrawText(to_string(i / step), xAxis + i + 2, 2);
	}

	for (int i = 0; i < 5000; i += gridWidth)
	{
		DrawText(to_string(-i / step), 2, yAxis - i + 2);
		DrawText(to_string(i / step), 2, yAxis + i + 2);
	}
}

void DrawPoints(vector<Point>& points)
{
	glColor3ub(255, 0, 0);
	glPointSize(10);
	glBegin(GL_POINTS);
	for (auto point : points)
		glVertex2f(xAxis + point.x * step, yAxis + point.y * step);
	glEnd();
	glPointSize(1);

	if (currentPoint != -1)
	{
		glColor3ub(0, 255, 128);
		glPointSize(10);
		glBegin(GL_POINTS);
		glVertex2f(xAxis + points[currentPoint].x * step, yAxis + points[currentPoint].y * step);
		glEnd();
		glPointSize(1);
	}
}

void DrawLines(vector<Point>& lines, GLubyte r, GLubyte g, GLubyte b)
{
	glLineWidth(2);
	glColor3ub(r, g, b);
	glBegin(GL_LINE_STRIP);
	for (auto line : lines)
		glVertex2f(xAxis + line.x * step, yAxis + line.y * step);
	glEnd();
	glLineWidth(1);
}

void Display()
{
	glClearColor(0.2, 0.2, 0.2, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH, GL_NICEST);

	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH, GL_NICEST);

	if (points.size() > 1)
	{
		if (s1)
		{
			DrawLines(points, 255, 255, 12);
		}
		if (s2)
		{
			LagrangePolynomial(points, spline1);
			DrawLines(spline1, 255, 153, 255);
		}
		if (s3)
		{
			LagrangeQuadratic(points, spline2);
			DrawLines(spline2, 17, 157, 220);
		}
		if (s4)
		{
			NaturalCubicSpline(points, spline3);
			DrawLines(spline3, 17, 255, 20);
		}
		if (s5)
		{
			BezierCurve(points, spline4);
			DrawLines(spline4, 156, 26, 84);
		}
	}
	else
		DrawText("Not enough points to draw spline", 2, Height - 17);

	if (Bpoints.size() > 0 && s6)
	{
		BSpline(Bpoints, spline5);

		if (Bpoints.size() < 4)
			DrawText("Not enough points to draw spline", 2, Height - 17);
		else
			DrawLines(spline5, 255, 255, 255);
	}

	DrawGrid();
	DrawLabels();
	DrawPoints(points);

	if (drawLabel && currentPoint != -1)
		DrawText(points[currentPoint].ToString(), xAxis + points[currentPoint].x * step, yAxis + points[currentPoint].y * step + 10);

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

void Mouse(int button, int state, int x, int y)
{
	int trueY = Height - y;

	if (state == GLUT_DOWN && button == GLUT_RIGHT_BUTTON)
	{
		points.push_back(Point((x - xAxis) / step, (trueY - yAxis) / step));
		Bpoints.push_back(points.back());
	}

	if (state == GLUT_DOWN && button == GLUT_LEFT_BUTTON)
	{
		if (points.size() > 0)
		{
			double gridX = (x - xAxis) / step;
			double gridY = (trueY - yAxis) / step;

			double min = (points[0].x - gridX) * (points[0].x - gridX) + (points[0].y - gridY) * (points[0].y - gridY);
			currentPoint = 0;
			for (int i = 1; i < points.size(); i++)
			{
				double dist = (points[i].x - gridX) * (points[i].x - gridX) + (points[i].y - gridY) * (points[i].y - gridY);

				if (dist < min)
				{
					currentPoint = i;
					min = dist;
				}
			}
		}
	}
}

void Keyboard(unsigned char key, int x, int y)
{
	int BcurrentPoint = 0;
	if (currentPoint != -1)
	{
		for (int i = 0; i < Bpoints.size(); i++)
			if ((Bpoints[i].y - points[currentPoint].y) * (Bpoints[i].y - points[currentPoint].y) + (Bpoints[i].x - points[currentPoint].x) * (Bpoints[i].x - points[currentPoint].x) < 1.0e-10)
				BcurrentPoint = i;

		if (key == 'w')
		{
			points[currentPoint].y += 0.1;
			Bpoints[BcurrentPoint].y += 0.1;
		}

		if (key == 's')
		{
			points[currentPoint].y -= 0.1;
			Bpoints[BcurrentPoint].y -= 0.1;
		}

		if (key == 'a')
		{
			points[currentPoint].x -= 0.1;
			Bpoints[BcurrentPoint].x -= 0.1;
		}

		if (key == 'd')
		{
			points[currentPoint].x += 0.1;
			Bpoints[BcurrentPoint].x += 0.1;
		}

		if (key == 'z')
		{
			if (currentPoint > 0)
				currentPoint--;
			else
				currentPoint = points.size() - 1;
		}

		if (key == 'x')
		{
			if (currentPoint < points.size() - 1)
				currentPoint++;
			else
				currentPoint = 0;
		}
	}

	if (key == 'c')
	{
		points.clear();
		Bpoints.clear();

		spline1.clear();
		spline2.clear();
		spline3.clear();
		spline4.clear();
		spline5.clear();

		currentPoint = -1;
	}


	if (key == '1')
		s1 = !s1;

	if (key == '2')
		s2 = !s2;

	if (key == '3')
		s3 = !s3;

	if (key == '4')
		s4 = !s4;

	if (key == '5')
		s5 = !s5;

	if (key == '6')
		s6 = !s6;

	if (key == 'i')
		yAxis -= 25;

	if (key == 'k')
		yAxis += 25;

	if (key == 'j')
		xAxis += 25;

	if (key == 'l')
		xAxis -= 25;

	if (key == 'y' && step < 200)
		step += 10;

	if (key == 'h' && step > 10)
		step -= 10;

	if (key == 'v')
		drawLabel = !drawLabel;

	if (key == 'b' && currentPoint != -1)
	{
		Bpoints.erase(Bpoints.begin() + BcurrentPoint);

		for (int i = 0; i < Bpoints.size(); i++)
			if ((Bpoints[i].y - points[currentPoint].y) * (Bpoints[i].y - points[currentPoint].y) + (Bpoints[i].x - points[currentPoint].x) * (Bpoints[i].x - points[currentPoint].x) < 1.0e-10)
				Bpoints.erase(Bpoints.begin() + currentPoint);

		points.erase(points.begin() + currentPoint);
		currentPoint = -1;
	}

	glutPostRedisplay();
}

void main(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);
	glutInitWindowSize(Width, Height);
	glutCreateWindow("Splines");
	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutMouseFunc(Mouse);
	glutKeyboardFunc(Keyboard);
	glutMainLoop();
}
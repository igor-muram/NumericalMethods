#include "GL\freeglut.h"
#include <vector>
using namespace std;

int Width = 1280, Height = 720;
int xAxis = Width / 2, yAxis = Height / 2;
int step = 50;

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

void Display()
{
	glClearColor(0.5, 0.5, 0.5, 1);
	glClear(GL_COLOR_BUFFER_BIT);
	glEnable(GL_POINT_SMOOTH);

	DrawAxis();

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
		yAxis -= 10;

	if (key == 's')
		yAxis += 10;

	if (key == 'a')
		xAxis += 10;

	if (key == 'd')
		xAxis -= 10;

	glutPostRedisplay();
}

void Mouse(int button, int state, int x, int y)
{
	if (state != GLUT_DOWN) return;
	glutPostRedisplay();
}

void main(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);
	glutInitWindowSize(Width, Height);
	glutCreateWindow("Визуализация");
	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutKeyboardFunc(Keyboard);
	glutMouseFunc(Mouse);
	glutMainLoop();
}
#include "Window.h"
#include "AxisComponent.h"
#include "FunctionComponent.h"
#include "ParametricComponent.h"
#include "DifferenceComponent.h"

#include <fstream>
using namespace std;

struct point
{
	point(double x, double y) : x(x), y(y) {};
	double x, y;
};

vector<point> points;

void DrawRect1(SDL_Renderer* renderer, int x, int y, int size)
{
	SDL_Rect rect;
	rect.h = size;
	rect.w = size;
	rect.x = x - size / 2;
	rect.y = y - size / 2;

	SDL_RenderFillRect(renderer, &rect);
}

void DrawPoints(SDL_Renderer* renderer, AxisComponent* comp)
{
	int xAxis = comp->GetXAxis();
	int yAxis = comp->GetYAxis();
	int step = comp->GetStep();

	SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
	for (int i = 0; i < points.size() - 1; i++)
		SDL_RenderDrawLineF(renderer, xAxis + points[i].x * step, yAxis - points[i].y * step, xAxis + points[i + 1].x * step, yAxis - points[i + 1].y * step);

	SDL_SetRenderDrawColor(renderer, 0, 150, 0, 255);
	for (int i = 0; i < points.size(); i++)
		DrawRect1(renderer, xAxis + points[i].x * step, yAxis - points[i].y * step, 4);
}

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

Window::Window(string title, int x, int y, int width, int height)
{
	ReadPoints("C:/input/points.txt");

	SDL_Init(SDL_INIT_VIDEO);
	this->width = width;
	this->height = height;
	window = SDL_CreateWindow("An SDL2 window", x, y, width, height, SDL_WINDOW_OPENGL);
	renderer = SDL_CreateRenderer(window, 0, 0);
	isRunning = true;
}

void Window::Destroy()
{
	SDL_DestroyWindow(window);
	SDL_Quit();
}

void Window::HandleEvents()
{
	SDL_Event event;
	while (SDL_PollEvent(&event))
	{
		switch (event.type)
		{
		case SDL_QUIT:
			isRunning = false;
			break;
		}
	}
}

void Window::ProcessKeyboard()
{
	const Uint8* state = SDL_GetKeyboardState(nullptr);
	if (state[SDL_SCANCODE_ESCAPE])
		isRunning = false;

	if (state[SDL_SCANCODE_W])
		axisComp->AddYOffset(1.5);

	if (state[SDL_SCANCODE_S])
		axisComp->AddYOffset(-1.5);

	if (state[SDL_SCANCODE_A])
		axisComp->AddXOffset(1.5);

	if (state[SDL_SCANCODE_D])
		axisComp->AddXOffset(-1.5);

	if (state[SDL_SCANCODE_P])
		axisComp->AddStep(0.2);

	if (state[SDL_SCANCODE_M])
		axisComp->AddStep(-0.2);

	if (state[SDL_SCANCODE_ESCAPE])
		isRunning = false;
}

void Window::Render()
{
	SDL_SetRenderDrawColor(renderer, 255, 255, 255, SDL_ALPHA_OPAQUE);
	SDL_RenderClear(renderer);

	diffComp->Draw(renderer);
	axisComp->Draw(renderer);

	for (auto comp : functionComps)
		comp->Draw(renderer);

	for (auto comp : parametricComps)
		comp->Draw(renderer);

	DrawPoints(renderer, axisComp);
	SDL_RenderPresent(renderer);
}
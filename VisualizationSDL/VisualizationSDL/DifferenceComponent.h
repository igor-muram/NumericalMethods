#pragma once
#include "Window.h"
#include "AxisComponent.h"

#include "SDL/SDL.h"

#include <functional>

class DifferenceComponent
{
public:
	DifferenceComponent(Window* window, vector<function<double(double, double)>> F, double interval);
	double Norm(double x, double y);
	void Draw(SDL_Renderer* renderer);

private:
	vector<function<double(double, double)>> F;
	double interval;
	Window* window;
	AxisComponent* axis;

	int xAxis;
	int yAxis;
	int step;
};


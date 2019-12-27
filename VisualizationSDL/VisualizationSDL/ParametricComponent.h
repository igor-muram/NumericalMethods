#pragma once
#include "Window.h"
#include "AxisComponent.h"

#include "SDL/SDL.h"

#include <functional>

using namespace std;

class ParametricComponent
{
public:
	ParametricComponent(Window* window, double a, double b, double h, function<double(double)> xt, function<double(double)> yt);
	void Draw(SDL_Renderer* renderer);

private:
	Window* window;
	AxisComponent* axis;
	double a, b, h;
	function<double(double)> xt, yt;
};


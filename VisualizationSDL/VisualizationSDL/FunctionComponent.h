#pragma once
#include "Window.h"
#include "AxisComponent.h"

#include "SDL/SDL.h"

#include <functional>

using namespace std;

class FunctionComponent
{
public:
	FunctionComponent(Window* window, double a, double b, double h, function<double(double)> f);
	void Draw(SDL_Renderer* renderer);

private:
	Window* window;
	AxisComponent* axis;
	double a, b, h;
	function<double(double)> f;
};


#pragma once
#include "SDL/SDL.h"
#include "Window.h"

class AxisComponent
{
public:
	AxisComponent(Window* window, int step);
	void Draw(SDL_Renderer* renderer);

	void AddXOffset(double offset) { xAxis += offset; }
	void AddYOffset(double offset) { yAxis += offset; }
	void AddStep(double value)
	{
		if (value > 0 && step < 200)
			step += value;

		if (value < 0 && step > 10)
			step += value;
	}

	double GetXAxis() { return xAxis; }
	double GetYAxis() { return yAxis; }
	double GetStep() { return step; }

	void SetStep(double step) { this->step = step; }
	void SetXAxis(double xAxis) { this->xAxis = xAxis; }
	void SetYAxis(double yAxis) { this->yAxis = yAxis; }

private:
	Window* window;
	double xAxis;
	double yAxis;
	double step;
};


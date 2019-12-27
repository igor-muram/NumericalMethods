#include "AxisComponent.h"

AxisComponent::AxisComponent(Window* window, int step)
{
	this->window = window;
	xAxis = window->GetWidth() / 2;
	yAxis = window->GetHeight() / 2;
	this->step = step;
}

void AxisComponent::Draw(SDL_Renderer* renderer)
{
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
	SDL_RenderDrawLine(renderer, -1000, yAxis, 2000, yAxis);

	for (int i = step; i < 1000; i += step)
	{
		SDL_RenderDrawLine(renderer, xAxis - i, yAxis, xAxis - i, yAxis - 3);
		SDL_RenderDrawLine(renderer, xAxis + i, yAxis, xAxis + i, yAxis - 3);
	}

	SDL_RenderDrawLine(renderer, xAxis, -1000, xAxis, 1000);

	for (int i = step; i < 1000; i += step)
	{
		SDL_RenderDrawLine(renderer, xAxis, yAxis + i, xAxis + 3, yAxis + i);
		SDL_RenderDrawLine(renderer, xAxis, yAxis - i, xAxis + 3, yAxis - i);
	}
}

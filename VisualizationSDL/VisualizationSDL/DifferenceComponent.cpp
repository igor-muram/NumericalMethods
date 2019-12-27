#include "DifferenceComponent.h"

void DrawRect(SDL_Renderer* renderer, int x, int y, int size)
{
	SDL_Rect rect;
	rect.h = size;
	rect.w = size;
	rect.x = x - size / 2;
	rect.y = y - size / 2;

	SDL_RenderFillRect(renderer, &rect);
}

DifferenceComponent::DifferenceComponent(Window* window, vector<function<double(double, double)>> F, double interval)
{
	this->axis = window->GetAxisComponent();
	this->window = window;
	this->interval = interval;
	this->F = F;
}

double DifferenceComponent::Norm(double x, double y)
{
	double sum = 0;
	for (auto f : F)
		sum += f(x, y) * f(x, y);
	return sqrt(sum);
}

void DifferenceComponent::Draw(SDL_Renderer* renderer)
{
	xAxis = axis->GetXAxis();
	yAxis = axis->GetYAxis();
	step = axis->GetStep();

	for (int i = -300; i < 300; i++)
	{
		for (int j = -300; j < 300; j++)
		{
			double x = i * 0.02;
			double y = j * 0.02;

			double diff = Norm(x, y);

			if (diff < 0.5)
			{
				SDL_SetRenderDrawColor(renderer, 50, 50, 50, 255);
				DrawRect(renderer, xAxis + x * step, yAxis - y * step, 4);
			}

			if (diff < interval && diff > 0.5)
			{
				SDL_SetRenderDrawColor(renderer, 80, 80, 80, 255);
				DrawRect(renderer, xAxis + x * step, yAxis - y * step, 4);
			}

			for (int i = 1; i < 50; i++)
				if (diff > i* interval&& diff < (i + 1) * interval)
				{
					double color = 80 + i * 10;
					if (color > 255) color = 255;
					SDL_SetRenderDrawColor(renderer, color, color, color, 255);
					DrawRect(renderer, xAxis + x * step, yAxis - y * step, 4);
				}
		}
	}
}

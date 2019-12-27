#include "FunctionComponent.h"

FunctionComponent::FunctionComponent(Window* window, double a, double b, double h, function<double(double)> f)
{
	this->window = window;
	this->a = a;
	this->b = b;
	this->h = h;
	this->f = f;
	axis = window->GetAxisComponent();
}

void FunctionComponent::Draw(SDL_Renderer* renderer)
{
	SDL_SetRenderDrawColor(renderer, 255, 0, 0, SDL_ALPHA_OPAQUE);

	int n = (b - a) / h;

	int xAxis = axis->GetXAxis();
	int yAxis = axis->GetYAxis();
	int step = axis->GetStep();

	for (int i = 0; i < n - 1; i++)
	{
		double x0 = a + i * h;
		double x1 = a + (i + 1) * h;

		double f0 = -f(x0);
		double f1 = -f(x1);

		SDL_RenderDrawLineF(
			renderer,
			xAxis + x0 * step,
			yAxis + f0 * step,
			xAxis + x1 * step,
			yAxis + f1 * step
		);
	}
}
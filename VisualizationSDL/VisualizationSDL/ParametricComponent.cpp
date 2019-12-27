#include "ParametricComponent.h"

ParametricComponent::ParametricComponent(Window* window, double a, double b, double h, function<double(double)> xt, function<double(double)> yt)
{
	this->window = window;
	this->a = a;
	this->b = b;
	this->h = h;
	this->xt = xt;
	this->yt = yt;
	axis = window->GetAxisComponent();
}

void ParametricComponent::Draw(SDL_Renderer* renderer)
{
	SDL_SetRenderDrawColor(renderer, 255, 0, 0, SDL_ALPHA_OPAQUE);

	int n = (b - a) / h;

	int xAxis = axis->GetXAxis();
	int yAxis = axis->GetYAxis();
	int step = axis->GetStep();

	for (int i = 0; i < n - 1; i++)
	{
		double t0 = a + i * h;
		double t1 = a + (i + 1) * h;

		double x0 = xt(t0);
		double x1 = xt(t1);

		double y0 = yt(t0);
		double y1 = yt(t1);

		SDL_RenderDrawLineF(
			renderer,
			xAxis + x0 * step,
			yAxis + y0 * step,
			xAxis + x1 * step,
			yAxis + y1 * step
		);
	}
}
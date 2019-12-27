#include "Window.h"
#include "AxisComponent.h"
#include "FunctionComponent.h"
#include "ParametricComponent.h"
#include "DifferenceComponent.h"

vector<function<double(double, double)>> F =
{
	[](double x, double y) { return (x + 2) * (x + 2) + (y - 1) * (y - 1) - 4; },
	[](double x, double y) { return (x - 1) * (x - 1) + (y - 1) * (y - 1) - 4; },
	[](double x, double y) { return y - x + 3.5; }
};

int main(int argc, char* argv[])
{
	Window window("Visualization", 200, 200, 1280, 720);
	window.AddAxisComponent(new AxisComponent(&window, 100));
	window.AddFunctionComponent(new FunctionComponent(&window, -10, 10, 0.1, [](double t) { return t + 1; }));
	window.AddParametricComponent(new ParametricComponent(&window, 0, 10, 0.1, [](double t) { return -2 + 2 * cos(t); }, [](double t) { return -1 + 2 * sin(t); }));
	window.AddParametricComponent(new ParametricComponent(&window, 0, 10, 0.1, [](double t) { return 2 + 2 * cos(t); }, [](double t) { return -1 + 2 * sin(t); }));
	window.AddDiffComponent(new DifferenceComponent(&window, F, 2.0));

	while (window.IsRunning())
	{
		window.HandleEvents();
		window.ProcessKeyboard();
		window.Render();
	}

	window.Destroy();
	return 0;
}
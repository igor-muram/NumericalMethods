#pragma once

#include "SDL/SDL.h"
#include <string>
#include <vector>

using namespace std;

class AxisComponent;
class FunctionComponent;
class ParametricComponent;
class DifferenceComponent;

class Window
{
public:
	Window(string title, int x, int y, int width, int height);
	void Destroy();
	void HandleEvents();
	void ProcessKeyboard();
	void Render();

	void AddAxisComponent(AxisComponent* comp)
	{
		if (axisComp == nullptr)
			axisComp = comp;
	}

	void AddFunctionComponent(FunctionComponent* comp) { functionComps.push_back(comp); }
	void AddParametricComponent(ParametricComponent* comp) { parametricComps.push_back(comp); }
	void AddDiffComponent(DifferenceComponent* comp) { diffComp = comp; }
	bool IsRunning() { return isRunning; }

	int GetWidth() { return width; }
	int GetHeight() { return height; }
	AxisComponent* GetAxisComponent() { return axisComp; }

private:
	int width = 0, height = 0;
	SDL_Window* window = nullptr;
	SDL_Renderer* renderer = nullptr;
	bool isRunning = false;
	AxisComponent* axisComp = nullptr;
	DifferenceComponent* diffComp = nullptr;

	vector<FunctionComponent*> functionComps;
	vector<ParametricComponent*> parametricComps;
};
#pragma once

struct Rect
{
	double x0 = 0.0, x1 = 0.0, z0 = 0.0, z1 = 0.0;
	Rect() = default;
	Rect(double x0, double x1, double z0, double z1) : x0(x0), x1(x1), z0(z0), z1(z1) {}
	Rect(const Rect& other) : x0(other.x0), x1(other.x1), z0(other.z0), z1(other.z1) {}
};

struct Element
{
	Rect rect;
	double value = 0.0;

	Element() : rect(Rect(0, 0, 0, 0)), value(0) {}
	Element(Rect rect, double value) : rect(rect), value(value) {}
};

struct Receiver
{
	double x = 0.0, z = 0.0;
	double value = 0.0;

	Receiver() = default;
	Receiver(double x, double z) : x(x), z(z) {}
	Receiver(double x, double z, double value) : x(x), z(z), value(value) {}
};
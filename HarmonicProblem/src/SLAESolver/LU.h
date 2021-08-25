#pragma once
#include "../stdafx.h"
#include "../Matrix.h"

namespace LU
{
	void LU(ProfileMatrix& A, std::vector<double>& x, std::vector<double>& b);
}
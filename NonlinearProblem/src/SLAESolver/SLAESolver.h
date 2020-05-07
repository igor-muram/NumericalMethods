#pragma once

#include <vector>

#include "../Matrix.h"

void LUDecomposition(Matrix& mat);
void Solve(Matrix& mat, std::vector<double>& x, std::vector<double>& b);
void Multiply(Matrix& mat, std::vector<double>& vec, std::vector<double>& res);

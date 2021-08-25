#pragma once

#include "pch.h"

#include <string>
#include <fstream>
#include <vector>

#include "Classes.h"

namespace CSV
{
	struct GridInfo
	{
		std::vector<Element>& elements;
		int ElementsCountX;
		int ElementsCountZ;
		double alpha;
		double difference;
	};

	void PrintElements(const std::string& filename, GridInfo info)
	{
		std::ofstream out(filename, std::ios::out);

		out << "Alpha;" << info.alpha << std::endl;
		out << "Difference: " << info.difference << std::endl;

		for (int i = 0; i < info.ElementsCountZ; i++)
		{
			for (int j = 0; j < info.ElementsCountX; j++)
			{
				out << info.elements[i * info.ElementsCountX + j].value << ";";
			}
			out << std::endl;
		}
		out.close();
	}

	void PrintElements(const std::string& filename, std::vector<Element>& elements, int ElementsCountX, int ElementsCountZ)
	{
		std::ofstream out(filename, std::ios::out);

		for (int i = 0; i < ElementsCountZ; i++)
		{
			for (int j = 0; j < ElementsCountX; j++)
			{
				out << elements[i * ElementsCountX + j].value << ";";
			}
			out << std::endl;
		}
		out.close();
	}

	void PrintElementsWithNumbers(const std::string& filename, std::vector<Element>& elements, int ElementsCountX, int ElementsCountZ)
	{
		std::ofstream out(filename, std::ios::out);

		for (int i = 0; i < ElementsCountZ; i++)
		{
			for (int j = 0; j < ElementsCountX; j++)
			{
				out << elements[i * ElementsCountX + j].value << " (" << i * ElementsCountX + j << ")" << ";";
			}
			out << std::endl;
		}
		out.close();
	}
}
#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <fstream>

class CSV
{
public:
	CSV(int isize, int jsize): isize(isize), jsize(jsize)
	{
		csv.resize(isize, std::vector<std::string>(jsize, ""));
	}

	template<typename T>
	void operator()(int i, int j, T value)
	{
		std::stringstream buffer;
		buffer << value;

		csv[i][j].append(buffer.str());
	}

	void Write(std::string filename)
	{
		std::ofstream out(filename);

		for (int i = 0; i < isize; i++)
		{
			for (int j = 0; j < jsize; j++)
			{
				out << csv[i][j] << ';';
			}

			out << std::endl;
		}

		out.close();
	}

private:

	int isize = 0, jsize = 0;
	std::vector<std::vector<std::string>> csv;
};
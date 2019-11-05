#pragma once

class Memory
{
public:
	Memory(int chunkSize, int numOfChunks) : 
		chunkSize(chunkSize),
		numOfChunks(numOfChunks),
		chunkNo(0),
		memory(new double[chunkSize * numOfChunks])
	{}

	void Allocate(double** ptr)
	{
		*ptr = &memory[chunkNo];
		chunkNo += chunkSize;
	}

	void Free(double** ptr)
	{
		*ptr = nullptr;
		chunkNo -= chunkSize;
	}

private:
	int chunkSize;
	int numOfChunks;
	int chunkNo;
	double* memory;
};
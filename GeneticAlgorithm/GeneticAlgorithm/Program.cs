using System;

namespace GeneticAlgorithm
{
	class Program
	{
		static void Main(string[] args)
		{

			double[] TrueGenes = new double[3] { 0.0, 1.0, 0.0 };

			GeneticAlgorithmInfo info = new GeneticAlgorithmInfo();
			info.PopulationCount = 1000;
			info.GenesCount = 3;
			info.TrueGenes = TrueGenes;
			info.MaxIterCount = 10000;
			info.Eps = 0.01;
			info.minGeneValue = -2.0;
			info.maxGeneValue = 2.0;
			info.MaxParentCount = 5;
			info.MutationProbability = 0.5;

			info.PointsCount = 1000;
			info.MinPoint = 0.0;
			info.MaxPoint = 2.0;

			GeneticAlgorithm algorithm = new GeneticAlgorithm(info);
			algorithm.Solve();

			Console.WriteLine();

			for (int i = 0; i < info.GenesCount; i++)
				Console.WriteLine($"{algorithm.BestGenotype[i]} - {TrueGenes[i]}");

			Console.WriteLine();
		}
	}

	
}

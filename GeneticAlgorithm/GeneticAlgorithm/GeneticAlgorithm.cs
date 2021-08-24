using System;
using System.Collections.Generic;

namespace GeneticAlgorithm
{
	struct GeneticAlgorithmInfo
	{ 
		public int PopulationCount { get; set; }
		public int GenesCount { get; set; }

		public double[] TrueGenes { get; set; }

		public int MaxIterCount { get; set; }
		public double Eps { get; set; }

		public double minGeneValue { get; set; }
		public double maxGeneValue { get; set; }

		public double MutationProbability { get; set; }
		public double MaxParentCount { get; set; }

		public int PointsCount { get; set; }
		public double MinPoint { get; set; }
		public double MaxPoint { get; set; }
	}

	class GeneticAlgorithm
	{
		public GeneticAlgorithmInfo Info { get; set; }

		Random Random = new Random();
		Population PreviousPopulation = null;
		Population CurrentPopulation = null;

		public List<double> BestFuctionalValues { get; set; } = new List<double>();
		public int Iteration { get; set; } = 0;

		public double[] BestGenotype { get; set; } = null;

		public GeneticAlgorithm(GeneticAlgorithmInfo info) => Info = info;

		public void Solve()
		{
			PreviousPopulation = CreateInitialPopulation();
			PreviousPopulation.SortIndividuals();
			BestFuctionalValues.Add(PreviousPopulation.Individuals[0].FunctionalValue);

			double bestFunctionalValue = BestFuctionalValues[0];
			Iteration = 0;
			while (Iteration < Info.MaxIterCount && bestFunctionalValue >= Info.Eps)
			{
				CurrentPopulation = CreateNewPopulation();
				CurrentPopulation.SortIndividuals();
				BestFuctionalValues.Add(CurrentPopulation.Individuals[0].FunctionalValue);
				bestFunctionalValue = CurrentPopulation.Individuals[0].FunctionalValue;

				PreviousPopulation = Selection();
				CurrentPopulation = null;
				Iteration++;

				Console.WriteLine(bestFunctionalValue);
			}

			BestGenotype = PreviousPopulation.Individuals[0].Genes;
		}

		Population CreateInitialPopulation()
		{	
			Population population = new Population(Info.PopulationCount);

			for (int i = 0; i < Info.PopulationCount; i++)
			{
				Individual individual = new Individual(Info.GenesCount);

				for (int j = 0; j < Info.GenesCount; j++)
					individual.Genes[j] = Info.minGeneValue + Random.NextDouble() * (Info.maxGeneValue - Info.minGeneValue);

				SetFunctionalValue(individual);
				population.Individuals[i] = individual;
			}

			return population;
		}

		Population CreateNewPopulation()
		{
			Population population = new Population(Info.PopulationCount);

			int n = 0;
			for (int i = 0; i < Info.MaxParentCount; i++)
			{
				for (int j = 0; j < Info.PopulationCount / Info.MaxParentCount; j++)
				{
					Individual father = PreviousPopulation.Individuals[i];
					Individual mother = PreviousPopulation.GetRandomIndividual(Random);

					Individual child = CreateChild(father, mother);
					Mutation(child);
					SetFunctionalValue(child);
					population.Individuals[n] = child;
					n++;
				}
			}

			return population;
		}

		Individual CreateChild(Individual mother, Individual father)
		{
			Individual child = new Individual(Info.GenesCount);

			int IndCross = Random.Next(Info.GenesCount);
			for (int i = 0; i < IndCross; i++)
				child.Genes[i] = father.Genes[i];

			for (int i = IndCross; i < Info.GenesCount; i++)
				child.Genes[i] = mother.Genes[i];

			return child;		
		}

		void Mutation(Individual individual)
		{
			double p = Random.NextDouble();
			if (p < Info.MutationProbability)
			{
				int i = Random.Next(Info.GenesCount);
				individual.Genes[i] = Info.minGeneValue + Random.NextDouble() * (Info.maxGeneValue - Info.minGeneValue);
			}
		}

		Population Selection()
		{
			Population population = new Population(Info.PopulationCount);

			int k = 0;
			int m = 0;
			for (int i = 0; i < Info.PopulationCount; i++)
			{
				double F1 = PreviousPopulation.Individuals[k].FunctionalValue;
				double F2 = CurrentPopulation.Individuals[m].FunctionalValue;
				if (F1 < F2)
				{
					population.Individuals[i] = PreviousPopulation.Individuals[k];
					k++;
				}
				else
				{
					population.Individuals[i] = CurrentPopulation.Individuals[m];
					m++;
				}
			}

			return population;
		}

		void SetFunctionalValue(Individual individual)
		{
			double result = 0.0;

			double h = (Info.MaxPoint - Info.MinPoint) / Info.PointsCount;
			for (int i = 0; i < Info.PointsCount; i++)
			{
				double point = Info.MinPoint + h * i;
				double value1 = Polynom(point, Info.TrueGenes);
				double value2 = Polynom(point, individual.Genes);
				result += (value1 - value2) * (value1 - value2) / (value1 * value1 + 1);
			}

			individual.FunctionalValue = Math.Sqrt(result);
		}

		double Polynom(double point, double[] coeffs)
		{
			double result = 0.0;
			for (int i = 0; i < Info.GenesCount; i++)
				result = result * point + coeffs[i];

			return result;
		}
	}

	class Individual
	{
		public double[] Genes { get; set; }
		public double FunctionalValue { get; set; }
		public Individual(int genesCount) => Genes = new double[genesCount];
	}

	class Population
	{
		public Individual[] Individuals = null;

		public Population(int populationCount) => Individuals = new Individual[populationCount];

		public void SortIndividuals() => Array.Sort(Individuals, (Individual x, Individual y) => x.FunctionalValue.CompareTo(y.FunctionalValue));

		public Individual GetRandomIndividual(Random random) => Individuals[random.Next(Individuals.Length)];
	}
}

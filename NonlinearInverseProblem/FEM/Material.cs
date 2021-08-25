namespace FEM
{
	public class Material
	{
		public double Sigma { get; set; } = 0.0;
		public double Gamma { get; set; } = 0.0;

		public Material(double sigma = 0.0, double gamma = 0.0)
		{
			Sigma = sigma;
			Gamma = gamma;
		}
	}
}

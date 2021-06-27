using System;

using BasicBFB.Testing;

namespace BasicBFB
{
	class Program
	{
		static void Main(string[] args)
		{
			// TODO: Upload to Github server

			GasifierParams param = ExampleData.gasifierParams();
			Pyrolysis pyro = new Pyrolysis(param);
			pyro.pyrolize();

			Console.WriteLine("\nWell I made it past the calls for pyrolysis...\n");

			double[] wTarCHO = pyro.tarCHO.w;
			double[] wGasCHO = pyro.dryGasCHO.w;
			double[] xGas = pyro.dryGasOut.x;

			string s = "";
			s = String.Format("Tar wC: {0:g3}  wH: {1:g3}  wO: {2:g3}",
								wTarCHO[0], wTarCHO[1], wTarCHO[2]);
			s += "\n";

			s += String.Format("Gas wC: {0:g3}  wH: {1:g3}  wO: {2:g3}",
								wGasCHO[0], wGasCHO[1], wGasCHO[2]);
			s += "\n\n";

			s += String.Format("Gas xCO:  {0:g3}   xCO2: {1:g3}", xGas[0], xGas[1]);
			s += "\n";
			s += String.Format("    xH2:  {0:g3}   xH2O: {1:g3}", xGas[2], xGas[3]);
			s += "\n\n";
			Console.WriteLine("s");
			Console.ReadKey();
		}
	}
}

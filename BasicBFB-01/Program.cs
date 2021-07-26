using System;
using System.Collections.Generic;

using BasicBFB.Model;
using BasicBFB.Testing;
using Common = BasicBFB.Model.Common;

namespace BasicBFB
{
	class Program
	{
		static void Main(string[] args)
		{
			GasifierParams param = ExampleData.gasifierParams();
			ReactorBed bed = new ReactorBed(ref param);
			Effluent bedEffluent = bed.calcEffluent();

			ReactorFreeboard freeboard = new ReactorFreeboard(param, bed.pyro, bedEffluent);
			Effluent finalEffluent = freeboard.calcEffluent();

			double pctGas = 100.0 * finalEffluent.productYields[0];
			double pctTar = 100.0 * finalEffluent.productYields[1];
			double pctChar = 100.0 * finalEffluent.productYields[2];

			double gasToFd = finalEffluent.yieldsVsDryBM[0];
			double tarToFd = finalEffluent.yieldsVsDryBM[1];
			double charToFd = finalEffluent.yieldsVsDryBM[2];

			double pctCO = 100.0 * finalEffluent.productGas.x[(int)Common.Component.CO];
			double pctCO2 = 100.0 * finalEffluent.productGas.x[(int)Common.Component.CO2];
			double pctCH4 = 100.0 * finalEffluent.productGas.x[(int)Common.Component.CH4];
			double pctH2 = 100.0 * finalEffluent.productGas.x[(int)Common.Component.H2];

			string dash80 = new String('=', 80);
			string report = "\n" + dash80 + "\n\n";

			report += "Full reactor model complete using fixed alpha = 0.9\n\n";

			report += "Product yields (dry mass basis):\n";
			report += $"\t--->  Gas:  {pctGas:F2} %\t\t{gasToFd:F3} kg/kg vs feed\n";
			report += $"\t--->  Tar:  {pctTar:F2} %\t\t{tarToFd:F3} kg/kg vs feed\n";
			report += $"\t--->  Char: {pctChar:F2} %\t\t{charToFd:F3} kg/kg vs feed\n\n";

			report += "Product gas composition (volume basis):\n";
			report += $"\t--->  CO:  {pctCO:F2} %\n";
			report += $"\t--->  CO2: {pctCO2:F2} %\n";
			report += $"\t--->  CH4: {pctCH4:F2} %\n";
			report += $"\t--->  H2:  {pctH2:F2} %\n";

			report += "\n" + dash80 + "\n\n";
			report += "Thank you and have a nice day!";

			Console.WriteLine(report);
			Console.ReadKey();
		}
	}
}

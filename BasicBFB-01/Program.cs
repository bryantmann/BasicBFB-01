using System;

using BasicBFB.Testing;

namespace BasicBFB
{
	class Program
	{
		static void Main(string[] args)
		{
			GasifierParams p = ExampleData.gasifierParams();
			CrossCheck tester = new CrossCheck(p);

			tester.calcGasFracs();
			tester.printComparisons();

			Console.WriteLine("\n\t\t\t\t- FIN -");
			Console.ReadKey();
		}
	}
}

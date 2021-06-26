using System;

using BasicBFB.Testing;

namespace BasicBFB
{
	class Program
	{
		static void Main(string[] args)
		{
			// TODO: Upload to Github server ASAP!!

			GasifierParams p = ExampleData.gasifierParams();
			CrossCheck tester = new CrossCheck(p);

			tester.calcGasFracs();

			// Arbitraty path for file
			string filepath = @"C:\Users\Bryan\temp\gasData.csv";
			string csvData = tester.dataToCSV();

			Utilities.saveTextAs(csvData, filepath);

			Console.WriteLine("File saved to bryan's temp directory as gasData.csv");
			Console.ReadKey();
		}
	}
}

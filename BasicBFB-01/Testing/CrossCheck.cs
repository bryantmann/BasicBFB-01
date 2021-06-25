using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;

namespace BasicBFB.Testing
{
	public class CrossCheck
	{
		public GasifierParams param { get; set; }
		public Pyrolysis pyro { get; set; }

		// Recall indices for C = 0, H = 1, and O = 2
		public double[,] gasFromH { get; set; }
		public double[,] gasFromC { get; set; }

		// 3D array to define constraints on composition
		//	- Outer rank - Holds the set of 2D arrays
		//  - 1st inner rank counts rows listing wC, wH and wO triplets
		//  - 2nd inner rank (columns) lists calculated values from given xH value
		public double[][,] constraints { get; set; }

		const int numComp = 3;
		const double wHmax = 0.25;
		public int numPoints { get; set; } = 100;


		// --------------------------------------------------------------------------
		// -------------------------- Constructors, etc -----------------------------
		// --------------------------------------------------------------------------

		public CrossCheck(GasifierParams param)
		{
			this.param = param;
			this.pyro = new Pyrolysis(param);
			this.gasFromH = new double[numPoints, numComp];
			this.gasFromC = new double[numPoints, numComp];

			this.constraints = new double[2][,];
			this.constraints[0] = new double[numPoints, numComp];
			this.constraints[1] = new double[numPoints, numComp];

			// Fill in values of wH that span range from 0 to 0.25
			const double deltaH = wHmax / ((double)numComp);
			double wHtemp = 0.0;
			for (int i = 0; i < numPoints; i++)
			{
				wHtemp += deltaH;
				gasFromH[i, 0] = 0.0;       // Carbon to be calculated
				gasFromH[i, 1] = wHtemp;
				gasFromH[i, 2] = 0.0;       // Oxygen, same deal

				constraints[0][i, 1] = wHtemp;
				constraints[1][i, 1] = wHtemp;
			}
		}

		public void calcGasFracs()
		{
			for (int i = 0; i < numPoints; i++)
			{
				gasFromH[i, 0] = pyro.pyroGasC(gasFromH[i, 1]);
				gasFromH[i, 2] = 1.0 - gasFromH[i, 0] - gasFromH[i, 1];

				// Use same values of carbon for cross-checks
				gasFromC[i, 0] = gasFromH[i, 0];
				gasFromC[i, 1] = pyro.pyroGasH(gasFromC[i, 0]);
				gasFromC[i, 2] = 1.0 - gasFromC[i, 0] - gasFromC[i, 1];
			}
		}


		// Index [0][j, k] contains low values, and [1][j, k] high values
		// After cross-check OK, set from one of the two sets of gas yields
		// Currently using gasFromH but we'll see
		public void setConstraints()
		{
			for (int i = 0; i < numPoints; i++)
			{
				double wH = constraints[0][i, 1];

				// Lower bounds (xC @ j = 0)
				constraints[0][i, 0] = 0.2727 + 1.9091 * wH;
				constraints[0][i, 2] = 0.5714 - 2.2857 * wH;

				// Upper bounds
				constraints[1][i, 0] = 0.4286 + 1.9091 * wH;
				constraints[1][i, 2] = 0.7273 - 2.2857 * wH;
			}
		}


		public List<double[]> validatedSets(double[,] wRaw)
		{
			List<double[]> wGood = new List<double[]>();

			for (int i = 0; i < numPoints; i++)
			{
				var rngC = (min: constraints[0][i, 0], max: constraints[1][i, 0]);
				var rngO = (min: constraints[0][i, 2], max: constraints[1][i, 2]);

				// The array below stores cross-check results for C and O
				// Indexed different - C is at 0 while O is at 1 (usually 2)
				bool[] xcheck = new bool[2];

				xcheck[0] = (wRaw[i, 0] > rngC.min) && (wRaw[i, 0] < rngC.max);
				xcheck[1] = (wRaw[i, 2] > rngO.min) && (wRaw[i, 2] < rngO.max);

				// If both xcheck values are true, then add to wGood
				if (xcheck[0] && xcheck[1])
				{
					double[] w = new double[numComp];
					w[0] = wRaw[i, 0];
					w[1] = wRaw[i, 1];
					w[2] = wRaw[i, 2];
					wGood.Add(w);
				}
			}

			return wGood;
		}


		public double[,] deltas(double[,] arr1, double[,] arr2)
		{
			int numRows = arr1.GetLength(0);
			int numCols = arr1.GetLength(1);
			double[,] d = new double[numRows, numCols];

			for (int i = 0; i < numRows; i++)
			{
				for (int j = 0; j < numCols; j++)
				{
					d[i, j] = Math.Abs(arr1[i, j] - arr2[i, j]);
				}
			}

			return d;
		}


		// Assume both arr1 and arr2 have same dimensions
		public double[] eqnRMSD(double[,] arr1, double[,] arr2)
		{
			int numRows = arr1.GetLength(0);
			int numCols = arr1.GetLength(1);
			double[] rmsd = new double[3] { 0.0, 0.0, 0.0 };
			double[] ssr = new double[3] { 0.0, 0.0, 0.0 };

			for (int i = 0; i < numRows; i++)
			{
				for (int j = 0; j < numCols; j++)
				{
					double r = arr1[i, j] - arr2[i, j];
					r = r * r;
					ssr[j] += r;
				}
			}

			// convert ssr to rmsd
			for (int j = 0; j < numCols; j++)
			{
				rmsd[j] = Math.Sqrt(ssr[j] / numRows);
			}

			return rmsd;
		}


		public void printComparisons()
		{
			double[,] dx = deltas(gasFromH, gasFromC);
			double[] rmsd = eqnRMSD(gasFromH, gasFromC);
			double[] wH = new double[numPoints];

			const double deltaH = wHmax / ((double)numComp);
			wH[0] = deltaH;
			for (int i = 1; i < numPoints; i++)
			{
				wH[i] = wH[i - 1] + deltaH;
			}

			string line = new string('—', 80);

			Console.WriteLine("\n\t\tABS DIFFERENCES BY INDEX");
			Console.WriteLine(line + "\n");

			for (int i = 0; i < numPoints; i++)
			{
				string ln = "  " + i.ToString() + "\t\t";
				ln += wH[i].ToString() + "\t\t";
				ln += dx[i, 0].ToString() + "\t";
				ln += dx[i, 1].ToString() + "\t";
				ln += dx[i, 2].ToString();
				Console.WriteLine(ln);
			}

			string ln2 = "\nOverall RMSDs:\t\t";
			ln2 += "C: " + rmsd[0].ToString() + "  ";
			ln2 += "H: " + rmsd[1].ToString() + "  ";
			ln2 += "O: " + rmsd[2].ToString() + "\n";
			Console.WriteLine(ln2);
			Console.WriteLine(line);
			Console.WriteLine();
		}


		// Prints list of points that pass validation tests
		//void printValidation(List<double[]> wGood) 
		//{
		//	int n = wGood.Count;
		//	string line = new string('—', 80);

		//	Console.WriteLine('\n');
		//	Console.WriteLine(line);
		//	Console.WriteLine("           CONSISTENCY BETWEEN METHODS");



		//	Console.WriteLine(line);
		//	Console.WriteLine("\VALIDATION RESULTS\n");
		//	Console.WriteLine($"Number of Points: {n}  (out of {numPoints} tested)");

		//}

	}
}
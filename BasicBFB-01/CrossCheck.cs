using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicBFB_01
{
	public class CrossCheck
	{
		public GasifierParams param { get; set; }
		public Pyrolysis pyro { get; set; }

		// Recall indices for C = 0, H = 1, and O = 2
		public double[,] gasFromH { get; set; }
		public double[,] gasFromC { get; set; }

		// 3D array to define constraints on composition
		//	- First rank holds the min and max limits (not real points)
		//  - Second rank counts rows listing wC, wH and wO triplets
		//  - Third rank (columns) lists calculated values from given xH value
		public double[,,] constraints { get; set; }

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
			this.constraints = new double[2, numPoints, numComp];

			// Fill in values of wH that span range from 0 to 0.25
			const double deltaH = wHmax / ((double) numComp);
			double wHtemp = 0.0;
			for (int i = 0; i < numPoints; i++)
			{
				wHtemp += deltaH;
				gasFromH[i, 0] = 0.0;       // Carbon to be calculated
				gasFromH[i, 1] = wHtemp;
				gasFromH[i, 2] = 0.0;       // Oxygen, same deal

				constraints[0, i, 1] = wHtemp;
				constraints[1, i, 1] = wHtemp;
			}
		}


		// Calculates oxygen mass fraction from wH and wC via argument
		private void addOxygen(double[,] arr) 
		{
			for (int i = 0; i < numPoints; i++)
			{
				arr[i, 2] = 1.0 - arr[i, 0] - arr[i, 1];
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


		// Index [0, j, k] contains low values, and [1, j, k] high values
		// After cross-check OK, set from one of the two sets of gas yields
		// Currently using gasFromH but we'll see
		public void setConstraints()
		{
			for (int j = 0; j < numPoints; j++)
			{
				double wH = constraints[0, j, 1];

				// Lower bounds (xC @ j = 0)
				constraints[0, j, 0] = 0.2727 + 1.9091 * wH;
				constraints[0, j, 2] = 0.5714 - 2.2857 * wH;

				// Upper bounds
				constraints[1, j, 0] = 0.4286 + 1.9091 * wH;
				constraints[1, j, 2] = 0.7273 - 2.2857 * wH;
			}
		}





	}
}
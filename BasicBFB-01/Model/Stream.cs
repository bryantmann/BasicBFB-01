﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using BasicBFB.Model.Common;


namespace BasicBFB.Model
{
	/* Stream class
	 *		Represents a stream going in or out of a unit operation
	 *		Components considered here are molecular species (tar and char exceptions)
	 *		Basis can be either mass or molar - indicated by isMolar flag 
	 *		
	 *		Component fractions are stored in the x array with the following mappings:
	 *		  INDEX		  SPECIES
	 *			0			CO
	 *			1			CO2
	 *			2			CH4		(all hydrocarbons lumped in with CH4)
	 *			3			H2
	 *			4			H2O
	 *			5			N2		(may not be used much)
	 *			6			O2		(may not be used much)
	 *			7			S		(may not be used much)
	 *			8			Tar		
	 *			
	 *			[9			Char	(removed)]
	 */


	public class Stream
	{
		public static int numComp = 9;
		public bool isMolar { get; private set; }        // Default of false means mass basis by default
		public double[] x { get; set; }

		public double p { get; set; }				// bar
		public double T { get; set; }               // °C

		private double _flowrate = 0.0;				// Mass flow rate (kg/s)
		public double flowrate
		{
			get => _flowrate;
			set
			{
				if (value < 0)
				{
					throw new ArgumentOutOfRangeException(
						$"{nameof(value)} must be non-negative");
				}
				_flowrate = value;
			}
		}

		static public string[] componentNames = new string[9] { "CO", "CO2", "CH4", "H2",
																"H2O", "N2", "O2", "S species",
																"Tar"};

		// --------------------------------------------------------------------------
		// ---------------------------- CONSTRUCTORS --------------------------------
		// --------------------------------------------------------------------------
			

		// Default constructor.  temp in units of Celsius, and press in unit of bars
		public Stream(double temp = 25.0, double press = 1.01325, bool isMolar = false)
		{
			x = new double[] {0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0};  
			this.isMolar = isMolar;
			this.p = press;
			this.T = temp;
		}


		// Alternate constructor for initializing with a known composition for x
		// User can input newX array of any length up to value of Stream.numComp
		public Stream(double[] newX, double temp = 25.0, double press = 1.01325, bool isMolar = false)
		{
			this.isMolar = isMolar;
			this.p = press;
			this.T = temp;
			this.x = new double[10];
			
			int n = newX.Length;
			if (n > Stream.numComp)
			{
				n = Stream.numComp;
			}

			for (int i = 0; i < n; i++)
			{
				this.x[i] = newX[i];
			}
		}


		// --------------------------------------------------------------------------
		// -------------------------------- METHODS ---------------------------------
		// --------------------------------------------------------------------------




		public void setX(double[] newX)
		{
			int n = newX.Length;
			//double sum = 0.0;

			if (n > Stream.numComp)
			{
				n = Stream.numComp;
			}

			for (int i = 0; i < n; i++)
			{
				this.x[i] = newX[i];
				//sum += x[i];
			}

			if (x.Length > newX.Length)
			{
				// Backfill last elements with zeros
				for (int i = newX.Length; i < x.Length; i++)
				{
					this.x[i] = 0.0;
				}
			}

		}

		// TODO: Trim unused methods and properties, starting in Stream


		public void normalizeGasFractions()
		{
			double sumX = 0.0;
			int n = Stream.numComp;        // Excludes the last two components, Tar and Char
			if (n <= 0)
			{
				throw new ArgumentOutOfRangeException(
					$"Not enough components, only {n} components left after dropping tar and char");
			}

			for (int i = 0; i < n; i++)
			{
				sumX += x[i];
			}

			if (sumX <= 0.0)
			{
				throw new ArgumentOutOfRangeException(
					$"Sum of component fractions must be > 0 to normalize; instead it is {sumX}");
			}

			for (int i = 0; i < n; i++)
			{
				x[i] /= sumX;
			}

			// Zero out tar and char mass fractions
			x[n] = 0.0;
			x[n + 1] = 0.0;
		}

		// Returns an array of length n = numComp - 2
		public double[] gasesNormalized()
		{
			if (Stream.numComp < 3)
			{
				throw new ArgumentOutOfRangeException(
					$"Not enough components, only {Stream.numComp} components in Stream");
			}

			double[] y = new double[Stream.numComp];
			double sumX = 0.0;
			int n = Stream.numComp;

			for (int i = 0; i < n; i++)
			{
				sumX += x[i];
			}

			if (sumX <= 0.0)
			{
				throw new ArgumentOutOfRangeException(
					$"Sum of component fractions must be > 0 to normalize; instead it is {sumX}");
			}

			for (int i = 0; i < n; i++)
			{
				y[i] = x[i] / sumX;
			}

			return y;
		}


		// Normalize all components in x including tar, char so they sum to 1.0
		public void normalizeAllFractions()
		{
			double sumX = 0.0;
			for (int i = 0; i < Stream.numComp; i++)
			{
				sumX += x[i];
			}

			if (sumX <= 0.0)
			{
				throw new ArgumentOutOfRangeException(
					$"Sum of component fractions must be > 0 to normalize; instead it is {sumX}");
			}

			for (int i = 0; i < Stream.numComp; i++)
			{
				x[i] /= sumX;
			}

		}


		// Returns an array of length Stream.numComp containing normalized fractions
		public double[] allNormalized()		{
			if (Stream.numComp < 3)
			{
				throw new ArgumentOutOfRangeException(
					$"Not enough components, only {Stream.numComp} components in Stream");
			}

			double[] allX = new double[Stream.numComp];
			double sumX = 0.0;

			for (int i = 0; i < Stream.numComp; i++)
			{
				sumX += x[i];
			}

			if (sumX <= 0.0)
			{
				throw new ArgumentOutOfRangeException(
					$"Sum of component fractions must be > 0 to normalize; instead it is {sumX}");
			}

			for (int i = 0; i < Stream.numComp; i++)
			{
				allX[i] = x[i] / sumX;
			}

			return allX;
		}


		// Checks if the sum if X equals 1.0 within a tolerance NORM_TOL
		public bool isNormalized()
		{
			bool normal;
			double d = 1.0 - sumX();
			normal = Math.Abs(d) < Const.NORM_TOL;
			return normal;
		}

		public double sumX()
		{
			double sum = 0.0;
			foreach (double xi in x)
			{
				sum += xi;
			}
			return sum;
		}


		// Returns the object's stored properties arranged in a csv table for export
		public string dataToCSV(string name)
		{
			string csv = "";
			string label = name + "," + "," + "Stream object" + "\n";
			string headerRow = "Species," + "wt frac" + "\n";
			csv += label + "\n" + headerRow;

			// Tabulate contents of x
			for (int i = 0; i < numComp; i++)
			{
				csv += componentNames[i] + "," + x.ToString() + "\n";
			}

			return csv;
		}


		// --------------------------------------------------------------------------
		// ---------------------Computed Properties ---------------------------------
		// --------------------------------------------------------------------------


		// This returns the avg MW of the stream ignoring Tar and Char (assumes it is 0)
		public double avgMW
		{
			get
			{
				double avg = 0.0;
				int n = numComp;   

				// Copy the first numComp-2 values from x into local variable y and normalize
				// Normalization excludes char and tar since they are mostly not gas phase
				double[] y = gasesNormalized();

				if (isMolar)
				{
					// Avg MW is the weighted sum of MW_i using normalized mole fractions
					for (int i = 0; i < n; i++)
					{
						avg += y[i] * MW.all[i];
					}
				}
				else
				{
					// Calculate avg MW by summing w_i / MW_i for all i and taking inverse
					for (int i = 0; i < n; i++)
					{
						avg += y[i] / MW.all[i];
					}
					avg = 1.0 / avg;
				}

				return avg;
			}
		}


	}
}


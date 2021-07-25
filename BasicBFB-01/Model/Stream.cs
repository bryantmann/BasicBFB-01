using System;
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

		// --------------------------------------------------------------------------------
		//									CONSTRUCTORS
		// --------------------------------------------------------------------------------

		// Default constructor.  temp in units of Celsius, and press in unit of bars
		public Stream(double temp = 600.0, double press = 1.01325, bool isMolar = false)
		{
			x = new double[] { 0.25, 0.25, 0.25, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			this.isMolar = isMolar;
			this.p = press;
			this.T = temp;
		}


		// Alternate constructor for initializing with a known composition for x
		// User can input newX array of any length up to value of Stream.numComp
		public Stream(double[] newX, double temp = 600.0, double press = 1.01325, bool isMolar = false)
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


		// --------------------------------------------------------------------------------
		//								BUILDER METHODS
		// --------------------------------------------------------------------------------

		// Creates a new stream from mass flow rates.  Defaults to mass basis
		static public Stream CreateFromMassRates(double[] mDot, double T, double p, 
												 bool isMolar = false) 
		{
			Stream stream = new Stream(T, p, isMolar);
			stream.mDot = mDot;
			return stream;
		}

		// Creates new stream from molar flow rates.  Defaults to molar basis
		static public Stream CreateFromMoleRates(double[] nDot, double T, double p,
												 bool isMolar = true)
		{
			Stream stream = new Stream(T, p, isMolar);
			stream.mDot = nDot;
			return stream;
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


		// --------------------------------------------------------------------------------
		//								COMPUTED PROPERTIES
		// --------------------------------------------------------------------------------

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


		// Component mass rates
		public double[] mDot
		{
			get
			{
				double[] rates = new double[numComp];
				for (int j = 0; j < numComp; j++)
				{
					rates[j] = x[j] * _flowrate;
				}

				if (isMolar)
				{
					for (int j = 0; j < numComp; j++)
					{
						rates[j] = rates[j] * MW.all[j];
					}
				}
				return rates;
			}

			set
			{
				double[] rates = new double[numComp];
				double sumRates = 0.0;
				int n = value.Length;
				if (n > numComp)
				{
					n = numComp;
				}
				for (int j = 0; j < n; j++)
				{
					rates[j] = isMolar ? (value[j] / MW.all[j]) : value[j];
					sumRates += rates[j];
				}

				flowrate = sumRates;
				for (int j = 0; j < n; j++)
				{
					x[j] = (sumRates > 0.0) ? (rates[j] / sumRates) : 0.0;
				}
			}
		}


		// Component molar rates
		public double[] nDot
		{
			get
			{
				double[] rates = new double[numComp];
				for (int j = 0; j < numComp; j++)
				{
					rates[j] = x[j] * _flowrate;
				}

				if (!isMolar)
				{
					for (int j = 0; j < numComp; j++)
					{
						rates[j] = rates[j] / MW.all[j];
					}
				}
				return rates;
			}

			set
			{
				double[] rates = new double[numComp];
				double sumRates = 0.0;
				int n = value.Length;
				if (n > numComp)
				{
					n = numComp;
				}
				for (int j = 0; j < n; j++)
				{
					rates[j] = isMolar ? value[j] : (value[j] * MW.all[j]);
					sumRates += rates[j];
				}

				flowrate = sumRates;
				for (int j = 0; j < n; j++)
				{
					x[j] = (sumRates > 0.0) ? (rates[j] / sumRates) : 0.0;
				}
			}
		}


		// Total mass rate, kg/s
		public double mDotTotal
		{
			get
			{
				if (isMolar)
				{
					double[] massRates = mDot;
					double sumRates = 0.0;
					foreach (double mi in massRates)
					{
						sumRates += mi;
					}
					return sumRates;
				}
				else
				{
					return _flowrate;
				}
			}

			set
			{
				if (isMolar)
				{
					flowrate = value / avgMW;
				}
				else
				{
					flowrate = value;
				}
			}
		}


		// Total molar rate
		public double nDotTotal
		{
			get
			{
				if (isMolar)
				{
					return _flowrate;
				}
				else
				{
					double[] moleRates = nDot;
					double sumRates = 0.0;
					foreach (double ni in moleRates)
					{
						sumRates += ni;
					}
					return sumRates;
				}
			}

			set
			{
				if (isMolar)
				{
					flowrate = value ;
				}
				else
				{
					flowrate = value * avgMW;
				}
			}
		}

		// Ideal gas density, g/cm3
		public double rho
		{
			get
			{
				double rhoMolar = (p * 1.0e5) / (Const.Rgas * (T + 273.15));   // mol/m3
				return rhoMolar * avgMW * 1.0e-3;
			}
		}
	}
}


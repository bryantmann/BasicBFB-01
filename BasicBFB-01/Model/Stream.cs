using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicBFB_01.Model
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
	 *			7			H2S		(may not be used much)
	 *			8			Tar		(may not be used much)
	 *			9			Char	(may not be used much)
	 */


	class Stream
	{
		public static int numComp = 10;
		public bool isMolar { get; private set; }        // Default of false means mass basis by default

		// ----------------------------------------------------------------------
		

		public double[] x { get; set; }

		public double p { get; set; }				// Pascals
		public double T { get; set; }               // Kelvin

		private double _flowrate = 0.0;     // Mass or molar flow rate (kg/s or mol/s)
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

		// This returns the avg MW of the stream ignoring Tar and Char (assumes it is 0)
		public double avgMW
		{
			get 
			{
				double sumX = 0.0;
				double avg = 0.0;
				double n = numComp - 2;		// Excludes the last two components, Tar and Char
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




				return avg;
			}
		}
	
		

		// Default constructor.  temp in units of Kelvin, and press in unit of Pascal
		public Stream(double temp = 298.15, double press = 101325.0, bool isMolar = false)
		{
			this.isMolar = isMolar;
			this.p = press;
			this.T = temp;

		}


		//private double[] _x = new double[] {0.2, 0.2, 0.2, 0.2, 0.2,	
		//									0.0, 0.0, 0.0, 0.0, 0.0};       // Component mass or mole fractions 

		//public double[] x
		//{
		//	get => _x;
		//	set => _x = value;
		//}
	}
}

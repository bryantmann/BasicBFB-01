using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicBFB.Model.Common
{
	public static class Const
	{
		public const double pi = 3.1415926535897932384626433832795028841971693993751;
		public const double Rgas = 8.31446261815324;    // J / mol.K
		public const double g = 9.81;					// m/s2

		public const double NORM_TOL = 1.0e-7;			// Normalized enough
	}


	public static class MW
	{
		public const double C = 12.011;
		public const double H = 1.00794;
		public const double N = 14.00674;
		public const double O = 15.9994;
		public const double S = 32.065;

		public const double CO = MW.C + MW.O;
		public const double CO2 = MW.C + 2.0 * MW.O;
		public const double CH4 = MW.C + 4.0 * MW.H;
		public const double H2 = 2.0 * MW.H;
		public const double H2O = 2.0 * MW.H + MW.O;
		public const double O2 = 2.0 * MW.O;
		public const double N2 = 2.0 * MW.N;
		public const double H2S = MW.S + 2.0 * MW.H;
		public const double Tar = 6.0 * MW.C + 6.0 * MW.H + MW.O;	// As phenol
		public const double Char = MW.C;							// As carbon

		/*	molWts is an Immutable array with MW values at same index as described in Stream class
		 *	It is to be used for convenience and indexing purposes
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
		 */

		public static double[] all = new double[] { MW.CO, MW.CO2, MW.CH4, MW.H2, MW.H2O,
													   MW.N2, MW.O2, MW.S, MW.Tar };
	}

}

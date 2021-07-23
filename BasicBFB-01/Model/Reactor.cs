using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicBFB.Model
{
	public class Reactor
	{
		// --------------------------------------------------------------------------------
		//									PARAMETERS
		// --------------------------------------------------------------------------------
		public GasifierParams param;
		public Pyrolysis pyro;

		public List<double> z = new List<double>();		
		public List<double[]> mDot = new List<double[]>();  // j = CO, CO2, CH4, H2, H2O, N2, O2, S, Tar
		public double mChar;


		// --------------------------------------------------------------------------------
		//									CONSTRUCTORS
		// --------------------------------------------------------------------------------
		public Reactor(GasifierParams param, double mCharGuess)
		{
			this.param = param;
			this.mChar = mCharGuess;

			this.pyro = new Pyrolysis(param);
			pyro.pyrolize();
		}
	}
}

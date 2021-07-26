using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using BasicBFB.Model.Common;

namespace BasicBFB.Model
{
	public static class ReactorODEs
	{
		// This static class allows ODE methods to be used with
		// delegate* function pointers 

		// ODE for ReactorBed's system of reactions
		internal static double[] dmBed(double z, in double[] mz, ReactorZone rxr)
		{
			double[] dmdz = new double[Stream.numComp];
			double[] gasRates = rxr.gasifyRates(z, mz);
			double eps = rxr.epsilon(z, mz);
			double U = rxr.gasU(mz);

			dmdz[(int)Component.CO] = rxr.dmCO(mz, gasRates, U, eps);
			dmdz[(int)Component.CO2] = rxr.dmCO2(mz, gasRates, U, eps);
			dmdz[(int)Component.CH4] = rxr.dmCH4(mz, gasRates, U, eps);
			dmdz[(int)Component.H2] = rxr.dmH2(mz, gasRates, U, eps);
			dmdz[(int)Component.H2O] = rxr.dmH2O(mz, gasRates, U, eps);
			dmdz[(int)Component.Tar] = rxr.dmTar(mz, eps, U);
			dmdz[(int)Component.N2] = 0.0;
			dmdz[(int)Component.O2] = 0.0;
			dmdz[(int)Component.S] = 0.0;

			return dmdz;
		}

		// Function pointer for use with ODESolver
		//unsafe static delegate*<double, in double[], ReactorFreeboard, double[]> fbOdePtr;

		// ODE for ReactorFreeboard's system of reactions
		internal static double[] dmFreeboard(double z, in double[] mz, ReactorZone rxr)
		{
			double[] dmdz = new double[Stream.numComp + 1];
			double[] gasRates = rxr.gasifyRates(z, mz);
			double U = rxr.gasU(mz);
			double eps = 1.0;       // Not actually used, for inheritance purposes

			dmdz[(int)Component.CO] = rxr.dmCO(mz, gasRates, U, eps);
			dmdz[(int)Component.CO2] = rxr.dmCO2(mz, gasRates, U, eps);
			dmdz[(int)Component.CH4] = rxr.dmCH4(mz, gasRates, U, eps);
			dmdz[(int)Component.H2] = rxr.dmH2(mz, gasRates, U, eps);
			dmdz[(int)Component.H2O] = rxr.dmH2O(mz, gasRates, U, eps);
			dmdz[(int)Component.Tar] = rxr.dmTar(mz, U, eps);
			dmdz[(int)Component.N2] = 0.0;
			dmdz[(int)Component.O2] = 0.0;
			dmdz[(int)Component.S] = 0.0;
			dmdz[(int)Component.Char] = rxr.dmChar(mz, gasRates, U);
			return dmdz;
		}
	}
}

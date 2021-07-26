using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using BasicBFB.Model.Common;


namespace BasicBFB.Model
{
	// NOTE: Assumes plug flow of char particles, with same velocity as gas phase

	public class ReactorFreeboard : ReactorZone
	{
		private Effluent bedEffluentIn { get; set; }		// Results from ReactorBed calcs

		private double Lfb => (zSpan[1] - zSpan[0]);		// Length of freeboard section


		// --------------------------------------------------------------------------------
		//									CONSTRUCTOR
		// --------------------------------------------------------------------------------
		public ReactorFreeboard(in GasifierParams param, in Pyrolysis pyro, in Effluent bedEffluent)
		{
			this.param = param;
			this.pyro = pyro;		// To be copied from ReactorBed.pyro

			bedEffluentIn = bedEffluent;
			effluent = new Effluent();
			zSpan = new double[2] { bedEffluentIn.z, param.Ltot };

			int N = Stream.numComp + 1;
			this.mz0 = new double[N];
			double[] mDotGasIn = bedEffluent.total.mDot;
			for (int j = 0; j < mDotGasIn.Length; j++)
			{
				mz0[j] = mDotGasIn[j];
			}
			mz0[N-1] = bedEffluent.mDotCharOut;

			zList = new List<double>(200);
			mDotList = new List<double[]>(200);
		}


		// --------------------------------------------------------------------------------
		//							MAIN CALCULATION ALGORITHMS
		// --------------------------------------------------------------------------------

		unsafe override public void solve()
		{
			//zList.Clear();
			//mDotList.Clear();
			zSpan[0] = bedEffluentIn.z;
			zSpan[1] = param.Ltot;

			delegate*<double, in double[], ReactorZone, double[]> fbOdePtr = &ReactorODEs.dmFreeboard;
			OdeSolver odeSolver = new OdeSolver();
			odeSolver.rk45b(fbOdePtr, zSpan, mz0, this, ref zList, ref mDotList);

			// Populate effluent container with results at z = Ltot
			int iMax = zList.Count - 1;     // Last index of zList and mDotList
			double[] mDotOut = mDotList[iMax];
			int iChar = mDotOut.Length - 1;
			effluent = new Effluent(mDotOut, zList[iMax], param, pyro, mDotOut[iChar]);
		}


		// --------------------------------------------------------------------------------
		//							ODE FOR REACTION SYSTEM
		// --------------------------------------------------------------------------------
		protected override double[] dmDot(double z, in double[] mz)
		{
			double[] dmdz = new double[Stream.numComp + 1];
			double[] gasRates = gasifyRates(z, mz);
			double U = gasU(mz);
			double eps = 1.0;		// Not actually used, for inheritance purposes

			dmdz[(int)Component.CO] = dmCO(mz, gasRates, U, eps);
			dmdz[(int)Component.CO2] = dmCO2(mz, gasRates, U, eps);
			dmdz[(int)Component.CH4] = dmCH4(mz, gasRates, U, eps);
			dmdz[(int)Component.H2] = dmH2(mz, gasRates, U, eps);
			dmdz[(int)Component.H2O] = dmH2O(mz, gasRates, U, eps);
			dmdz[(int)Component.Tar] = dmTar(mz, U, eps);
			dmdz[(int)Component.N2] = 0.0;
			dmdz[(int)Component.O2] = 0.0;
			dmdz[(int)Component.S] = 0.0;
			dmdz[(int)Component.Char] = dmChar(mz, gasRates, U);
			return dmdz;
		}


		// --------------------------------------------------------------------------------
		//						REACTION RATE CALCS - HELPERS
		// --------------------------------------------------------------------------------

		// Gasification reaction rates returned in array
		// Order is rC1, rC2, rC3, rSMR, rWGS
		// Indexed the same as Reaction enumeration raw values
		internal override double[] gasifyRates(double z, in double[] mz)
		{
			double[] rates = new double[5];		// rC1, rC2, rC3, rSMR, rWGS
			double Tk = param.T + 273.15;       // T in Kelvin

			double pCO2 = partialP(Component.CO2, mz);
			double pH2O = partialP(Component.H2O, mz);
			double pH2  = partialP(Component.H2, mz);

			double[] cGas = new double[5];      // CO, CO2, CH4, H2, H2O
			cGas[0] = conc(Component.CO, mz);
			cGas[1] = conc(Component.CO2, mz);
			cGas[2] = conc(Component.CH4, mz);
			cGas[3]  = conc(Component.H2, mz);
			cGas[4] = conc(Component.H2O, mz);

			rates[0] = GasRxn.rateC1(pCO2, Tk);
			rates[1] = GasRxn.rateC2(pH2O, Tk);
			rates[2] = GasRxn.rateC3(pH2, Tk);
			rates[3] = GasRxn.rateSMR(cGas[2], cGas[4], Tk);
			rates[4] = GasRxn.rateWGS(cGas[0], cGas[1], cGas[3], cGas[4], Tk);

			return rates;
		}


		// CO rate expression (dmCO/dz)
		internal override double dmCO(in double[] mz, in double[] gasRates, double U, double eps)
		{
			double mDotChar = mz[(int)Component.Char];

			// Gasification term
			double gTerm1 = gasRates[(int)Reaction.SMR] - gasRates[(int)Reaction.WGS];
			double gTerm2 = mDotChar / (MW.Char * U * param.Axs);
			gTerm2 *= 2.0 * gasRates[(int)Reaction.C1] + gasRates[(int)Reaction.C2];
			double dmGas = (gTerm1 + gTerm2) * param.Axs * MW.CO;	

			// Tar cracking term
			double mzTar = mz[(int)Component.Tar];
			double dmTar = pyro.kPyros[3] * (1.0 - pyro.gammaTarInert); 
			dmTar *= (mzTar / U) * Arrhenius.wTarPyro[(int)Component.CO];

			return dmTar + dmGas;
		}


		// CO2 rate expression (dmCO2/dz)
		internal override double dmCO2(in double[] mz, in double[] gasRates, double U, double eps)
		{
			double mDotChar = mz[(int)Component.Char];

			// Gasification term
			double gTerm1 = gasRates[(int)Reaction.WGS];
			double gTerm2 = mDotChar / (MW.Char * U * param.Axs);
			gTerm2 *= -gasRates[(int)Reaction.C1];
			double dmGas = (gTerm1 + gTerm2) * param.Axs * MW.CO2;  // A major term

			// Tar cracking term
			double mzTar = mz[(int)Component.Tar];
			double dmTar = pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);    // A major term
			dmTar *= (mzTar / U) * Arrhenius.wTarPyro[(int)Component.CO2];

			return dmTar + dmGas;
		}


		// CH4 rate expression (dmCH4/dz)
		internal override double dmCH4(in double[] mz, in double[] gasRates, double U, double eps)
		{
			double mDotChar = mz[(int)Component.Char];

			// Gasification term
			double gTerm1 = -gasRates[(int)Reaction.SMR];
			double gTerm2 = mDotChar / (MW.Char * U * param.Axs);
			gTerm2 *= gasRates[(int)Reaction.C3];
			double dmGas = (gTerm1 + gTerm2) * param.Axs * MW.CH4;

			// Tar cracking term
			double mzTar = mz[(int)Component.Tar];
			double dmTar = pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);  
			dmTar *= (mzTar / U) * Arrhenius.wTarPyro[(int)Component.CH4];

			return dmTar + dmGas;
		}

		// H2 rate expression (dmH2/dz)
		internal override double dmH2(in double[] mz, in double[] gasRates, double U, double eps)
		{
			double mDotChar = mz[(int)Component.Char];

			// Gasification term
			double gTerm1 = 3.0 * gasRates[(int)Reaction.SMR] + gasRates[(int)Reaction.WGS];
			double gTerm2 = mDotChar / (MW.Char * U * param.Axs);
			gTerm2 *= gasRates[(int)Reaction.C2] - 2.0 * gasRates[(int)Reaction.C3];
			double dmGas = (gTerm1 + gTerm2) * param.Axs * MW.H2;  // A major term

			// Tar cracking term
			double mzTar = mz[(int)Component.Tar];
			double dmTar = pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);    // A major term
			dmTar *= (mzTar / U) * Arrhenius.wTarPyro[(int)Component.H2];

			return dmTar + dmGas;
		}

		// H2O rate expression (dmH2O/dz)
		internal override double dmH2O(in double[] mz, in double[] gasRates, double U, double eps)
		{
			double mDotChar = mz[(int)Component.Char];

			// Gasification term
			double gTerm1 = -(gasRates[(int)Reaction.SMR] + gasRates[(int)Reaction.WGS]);
			double gTerm2 = mDotChar / (MW.Char * U * param.Axs); ;
			gTerm2 *= -gasRates[(int)Reaction.C2];
			double dmGas = (gTerm1 + gTerm2) * param.Axs * MW.H2O;

			return dmGas;
		}


		// Tar rate expression (dmTar/dz)
		internal override double dmTar(in double[] mz, double U, double eps)
		{
			double dmTar = -pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);
			dmTar *= mz[(int)Component.Tar] / U;
			return dmTar;
		}


		// Char rate expression (dmChar/dz)
		internal override double dmChar(in double[] mz, double[] gasRates, double U)
		{
			double mDotChar = mz[(int)Component.Char];
			double dmChar = -gasRates[(int)Reaction.C1] - gasRates[(int)Reaction.C2] 
							- gasRates[(int)Reaction.C3];
			dmChar *= mDotChar / U;

			return dmChar;
		}
	}
}

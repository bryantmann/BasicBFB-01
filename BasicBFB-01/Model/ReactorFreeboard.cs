using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using BasicBFB.Model.Common;


namespace BasicBFB.Model
{
	// HACK: Assumes plug flow of char particles, with same velocity as gas phase
	public class ReactorFreeboard
	{
		public GasifierParams param;
		public Pyrolysis pyro;
		public Effluent bedEffluentIn { get; init; }
		
		// Final reactor effluent
		public Effluent effluent { get; private set; }

		private double[] zSpan { get; init; }     // {Lbed, Ltot}
		private double[] mz0 { get; init; }       // Initial condition = Mass rates from reactor bed out

		public List<double> zList { get; private set; }
		public List<double[]> mDotList { get; private set; }
		// j = CO, CO2, CH4, H2, H2O, N2, O2, S, Tar, Char

		public double Lfb => (zSpan[1] - zSpan[0]);		// Length of freeboard section


		// --------------------------------------------------------------------------------
		//									CONSTRUCTOR
		// --------------------------------------------------------------------------------
		public ReactorFreeboard(GasifierParams param, Pyrolysis pyro, Effluent bedEffluent)
		{
			this.param = param;
			this.pyro = pyro;		// Copied from ReactorBed.pyro
			//this.pyro = new Pyrolysis(param);
			this.bedEffluentIn = bedEffluent;
			this.effluent = new Effluent();
			this.zSpan = new double[2] { bedEffluentIn.z, param.Ltot };

			this.mz0 = bedEffluent.total.mDot;
			this.zList = new List<double>(200);
			this.mDotList = new List<double[]>(200);
		}


		// --------------------------------------------------------------------------------
		//							MAIN CALCULATION ALGORITHMS
		// --------------------------------------------------------------------------------
		public Effluent calcEffluent()
		{
			solve();
			return effluent;
		}

		// TODO: Fill in this dude 
		public void solve()
		{
			// Code goes here
		}


		// --------------------------------------------------------------------------------
		//							ODE FOR REACTION SYSTEM
		// --------------------------------------------------------------------------------
		private double[] dmDot(double z, in double[] mz)
		{
			double[] dmdz = new double[Stream.numComp + 1];
			double[] gasRates = gasifyRates(z, mz);
			double U = gasU(mz);

			dmdz[(int)Component.CO] = dmCO(mz, gasRates, U);
			dmdz[(int)Component.CO2] = dmCO2(mz, gasRates, U);
			dmdz[(int)Component.CH4] = dmCH4(mz, gasRates, U);
			dmdz[(int)Component.H2] = dmH2(mz, gasRates, U);
			dmdz[(int)Component.H2O] = dmH2O(mz, gasRates, U);
			dmdz[(int)Component.Tar] = dmTar(mz, U);
			dmdz[(int)Component.N2] = 0.0;
			dmdz[(int)Component.O2] = 0.0;
			dmdz[(int)Component.S] = 0.0;
			dmdz[(int)Component.Char] = dmChar(mz, gasRates, U);
			return dmdz;
		}


		// --------------------------------------------------------------------------------
		//							  GENERAL HELPER FUNCTIONS
		// --------------------------------------------------------------------------------

		// Gas phase mole fraction for specified component
		// at z position corresponding to specified zIndex
		private double yMolar(Component comp, in double[] mz)
		{
			double sumMoles = 0.0;

			for (int j = 0; j < Stream.numComp; j++)
			{
				sumMoles += mz[j] / MW.all[j];
			}

			double y = mz[(int)comp] / MW.all[(int)comp];
			y /= sumMoles;

			return y;
		}

		// Partial pressure of specified component in gas phase 
		// at z position for specified index (bara)
		private double partialP(Component comp, in double[] mz)
		{
			return param.p * yMolar(comp, mz);
		}

		// Concentration of specified component in gas at zIndex [=] mol/m3
		private double conc(Component comp, in double[] mz)
		{
			double Tkelvin = param.T + 273.15;
			double c = 1.0e5 * partialP(comp, mz) / (Const.Rgas * Tkelvin);
			return c;
		}

		// Superficial gas phase velocity at z[zIndex]
		private double gasU(double[] mz)
		{
			double Tkelvin = param.T + 273.15;
			double u = 0.0;

			for (int j = 0; j < Stream.numComp; j++)
			{
				u += mz[j] / MW.all[j];
			}

			u *= (Const.Rgas * Tkelvin) / (1.0e5 * param.p * param.Axs);
			return u;
		}


		// --------------------------------------------------------------------------------
		//						REACTION RATE CALCS - HELPERS
		// --------------------------------------------------------------------------------

		// Gasification reaction rates returned in array
		// Order is rC1, rC2, rC3, rSMR, rWGS
		// Indexed the same as Reaction enumeration raw values
		private double[] gasifyRates(double z, in double[] mz)
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
		private double dmCO(in double[] mz, in double[] gasRates, double U)
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
		private double dmCO2(in double[] mz, in double[] gasRates, double U)
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
		private double dmCH4(in double[] mz, in double[] gasRates, double U)
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
		private double dmH2(in double[] mz, in double[] gasRates, double U)
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
		private double dmH2O(in double[] mz, in double[] gasRates, double U)
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
		private double dmTar(in double[] mz, double U)
		{
			double dmTar = -pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);
			dmTar *= mz[(int)Component.Tar] / U;
			return dmTar;
		}


		// Char rate expression (dmChar/dz)
		private double dmChar(in double[] mz, double[] gasRates, double U)
		{
			double mDotChar = mz[(int)Component.Char];
			double dmChar = -gasRates[(int)Reaction.C1] - gasRates[(int)Reaction.C2] 
							- gasRates[(int)Reaction.C3];
			dmChar *= mDotChar / U;

			return dmChar;
		}
	}
}

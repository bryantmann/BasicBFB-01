using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using BasicBFB.Model.Common;

namespace BasicBFB.Model
{
	// Well I fucked this one up.
	public class ReactorBed2
	{
		// --------------------------------------------------------------------------------
		//									PARAMETERS
		// --------------------------------------------------------------------------------
		public GasifierParams param;
		public Pyrolysis pyro;
		Stream pyroGas;

		public List<double> z = new List<double>();		
		public List<double[]> mDot = new List<double[]>();  // j = CO, CO2, CH4, H2, H2O, N2, O2, S, Tar

		public double mChar;
		public double Lbed;
		double epsAvg;

		const double TOL = 1.0e-3;		// Relative tolerance for mChar and Lbed convergence


		// --------------------------------------------------------------------------------
		//									CONSTRUCTOR
		// --------------------------------------------------------------------------------
		public ReactorBed2(GasifierParams param, double mCharGuess)
		{
			this.param = param;
			this.mChar = mCharGuess;

			this.pyro = new Pyrolysis(param);
			pyro.pyrolize();
			mixPyroGas();
			setInitialConditions();

			// Guess avgEps, then calculate Lbed using that
			epsAvg = 0.7;
			Lbed = param.L0 * (1.0 - param.epsMF) / (1.0 - epsAvg);
		}


		// --------------------------------------------------------------------------------
		//							MAIN CALCULATION ALGORITHM
		// --------------------------------------------------------------------------------
		public void solve()
		{
			double[] y0 = (double[]) mDot[0].Clone();
			double[] zSpan = { 0.0, Lbed };
			
		}



		// --------------------------------------------------------------------------------
		//							ODE FOR REACTION SYSTEM
		// --------------------------------------------------------------------------------
		private double[] dmDot(int zIndex)
		{
			double[] dmdz = new double[Stream.numComp];
			double[] gasRates = gasifyRates(zIndex);
			double eps = epsilon(zIndex);
			double U = gasU(zIndex);

			dmdz[(int)Component.CO] = dmCO(zIndex, gasRates, eps, U);
			dmdz[(int)Component.CO2] = dmCO2(zIndex, gasRates, eps, U);
			dmdz[(int)Component.CH4] = dmCH4(zIndex, gasRates, eps, U);
			dmdz[(int)Component.H2] = dmH2(zIndex, gasRates, eps, U);
			dmdz[(int)Component.H2O] = dmH2O(zIndex, gasRates, eps);
			dmdz[(int)Component.Tar] = dmTar(zIndex, eps, U);
			dmdz[(int)Component.N2] = 0.0;
			dmdz[(int)Component.O2] = 0.0;
			dmdz[(int)Component.S] = 0.0;

			return dmdz;
		}


		// --------------------------------------------------------------------------------
		//							CHAR HOLDUP CALCULATIONS
		// --------------------------------------------------------------------------------

		// Calculates new value for mChar and updates this object's mChar parameter value
		private void updateMChar()
		{
			double tauP = meanTauP();
			Del2 rCsum = sumCharRates;
			double integral = Solver.integrate(rCsum, z) / Lbed;
			double denom = (1.0 / tauP) + integral;
			mChar = pyro.charYield * pyro.dryFeedIn.flow / denom;
		}

		// d(tauP)/dz - derivative of particle residence time wrt z as a function of zIndex
		private double dtauP(int zIndex)
		{
			double eps = epsilon(zIndex);
			double U = gasU(zIndex);
			return (1.0 - eps) / U;
		}

		// Mean solids residence time tauP
		private double meanTauP()
		{
			Del2 dtau = dtauP;
			return Solver.integrate(dtau, z);
		}

		// Function used as integrand in mChar calculation method
		private double sumCharRates(int zIndex)
		{
			double Tk = param.T + 273.15;       // T in Kelvin

			double pCO2 = partialP(Component.CO2, zIndex);
			double pH2O = partialP(Component.H2O, zIndex);
			double pH2 = partialP(Component.H2, zIndex);

			double rC1 = GasRxn.rateC1(pCO2, Tk);
			double rC2 = GasRxn.rateC2(pH2O, Tk);
			double rC3 = GasRxn.rateC3(pH2, Tk);
			return rC1 + rC2 + rC3;
		}


		// --------------------------------------------------------------------------------
		//								SETUP HELPER FUNCTIONS
		// --------------------------------------------------------------------------------

		// Combine pyroGas = dryGasOut + badGasOut
		private void mixPyroGas()
		{
			this.pyroGas = new Stream(param.T, param.p);
			double[] f = new double[Stream.numComp];
			double sumf = 0.0;
			for (int i = 0; i < Stream.numComp; i++)
			{
				f[i] = pyro.dryGasOut.flowrate * pyro.dryGasOut.x[i]
						+ pyro.badGasOut.flowrate * pyro.dryGasOut.x[i];
				sumf += f[i];
			}

			pyroGas.flowrate = sumf;

			for (int i = 0; i < Stream.numComp; i++)
			{
				pyroGas.x[i] = f[i] / sumf;
			}
		}

		// Add first point to z = 0, and add initial conditions to mDot
		private void setInitialConditions()
		{
			z.Add(0.0);
			double[] mDot0 = new double[Stream.numComp];
			for (int i = 0; i < Stream.numComp; i++)
			{
				// Steam components are nonzero at z = 0, rest contribute nothing
				mDot0[i] = param.steamIn.flowrate * param.steamIn.x[i];
				mDot0[i] *= MW.all[i] / param.steamIn.avgMW;
			}
			mDot.Add(mDot0);
		}


		// --------------------------------------------------------------------------------
		//							  GENERAL HELPER FUNCTIONS
		// --------------------------------------------------------------------------------

		// Gas phase mole fraction for specified component
		// at z position corresponding to specified zIndex
		private double yMolar(Component comp, int zIndex)
		{
			double[] mz = this.mDot[zIndex];
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
		private double partialP(Component comp, int zIndex)
		{
			return param.p * yMolar(comp, zIndex);
		}

		// Concentration of specified component in gas at zIndex [=] mol/m3
		private double conc(Component comp, int zIndex)
		{
			double Tkelvin = param.T + 273.15;
			double c = 1.0e5 * partialP(comp, zIndex) / (Const.Rgas * Tkelvin);
			return c;
		}

		// Superficial gas phase velocity at z[zIndex]
		private double gasU(int zIndex)
		{
			double[] mz = this.mDot[zIndex];
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
		//						BUBBLING BED PROPERTY CALCULATIONS
		// --------------------------------------------------------------------------------

		// Estimated bubble size db as a function of height (per zIndex)
		private double bubbleSize(int zIndex)
		{
			double db = 0.0;
			double zi = z[zIndex];
			double deltaU = gasU(zIndex) - param.U0;
			double zTerm = zi + 4.0 * Math.Sqrt((param.Axs / ((double)param.Nor)));
			db = 0.54 * Math.Pow(Const.g, -0.2) * Math.Pow(deltaU, 0.4);
			db *= Math.Pow(zTerm, 0.8);
			return db;
		}

		// Local bed voidage eps(z) based on zIndex
		private double epsilon(int zIndex)
		{
			double deltaU = gasU(zIndex) - param.U0;
			double dbTerm = 0.711 * Math.Sqrt(Const.g * bubbleSize(zIndex));
			double denom = 1.0 + deltaU / dbTerm;
			return 1.0 - (1.0 - param.epsMF) / denom;
		}

		// Average voidage of the whole fluidized bed
		private double avgEpsilon()
		{
			Del2 eps = epsilon;
			double epsIntegral = Solver.integrate(eps, z);
			return epsIntegral / Lbed;
		}

		// Calculate height of bubbling bed and set value of Lbed
		private void updateBedHeight()
		{
			epsAvg = avgEpsilon();
			Lbed = param.L0 * (1.0 - param.epsMF) / (1.0 - epsAvg);
		}


		// --------------------------------------------------------------------------------
		//						REACTION RATE CALCS - HELPERS
		// --------------------------------------------------------------------------------

		// Gasification reaction rates returned in array
		// Order is rC1, rC2, rC3, rSMR, rWGS
		// Indexed the same as Reaction enumeration raw values
		private double[] gasifyRates(int zIndex)
		{
			double[] rates = new double[5];		// rC1, rC2, rC3, rSMR, rWGS
			double Tk = param.T + 273.15;       // T in Kelvin

			double pCO2 = partialP(Component.CO2, zIndex);
			double pH2O = partialP(Component.H2O, zIndex);
			double pH2  = partialP(Component.H2, zIndex);

			double[] cGas = new double[5];      // CO, CO2, CH4, H2, H2O
			cGas[0] = conc(Component.CO, zIndex);
			cGas[1] = conc(Component.CO2, zIndex);
			cGas[2] = conc(Component.CH4, zIndex);
			cGas[3]  = conc(Component.H2, zIndex);
			cGas[4] = conc(Component.H2O, zIndex);

			rates[0] = GasRxn.rateC1(pCO2, Tk);
			rates[1] = GasRxn.rateC2(pH2O, Tk);
			rates[2] = GasRxn.rateC3(pH2, Tk);
			rates[3] = GasRxn.rateSMR(cGas[2], cGas[4], Tk);
			rates[4] = GasRxn.rateWGS(cGas[0], cGas[1], cGas[3], cGas[4], Tk);

			return rates;
		}


		// CO rate expression (dmCO/dz)
		private double dmCO(int zIndex, in double[] gasRates, double eps, double U)
		{
			double gTerm1 = eps * (gasRates[(int)Reaction.SMR] - gasRates[(int)Reaction.WGS]);
			
			double gTerm2 = (1.0 - eps) * (mChar / MW.Char);
			gTerm2 /= param.Axs * Lbed * (1.0 - eps);
			gTerm2 *= 2.0 * gasRates[(int)Reaction.C1] + gasRates[(int)Reaction.C2];
			
			double dmGas = (gTerm1 + gTerm2) * param.Axs * MW.CO;	// A major term

			double term1b = pyro.gasYield * pyro.dryFeedIn.flow;    // A major term
			term1b *= pyro.dryGasOut.x[(int)Component.CO] / Lbed;

			double[] mz = mDot[zIndex];
			double mzTar = mz[(int)Component.Tar];

			double term2b = pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);    // A major term
			term2b *= (mzTar * eps / U) * Arrhenius.wTarPyro[(int)Component.CO];

			return term1b + term2b + dmGas;
		}


		// CO2 rate expression (dmCO2/dz)
		private double dmCO2(int zIndex, in double[] gasRates, double eps, double U)
		{
			double gTerm1 = eps * gasRates[(int)Reaction.WGS];

			double gTerm2 = (1.0 - eps) * (mChar / MW.Char);
			gTerm2 /= param.Axs * Lbed * (1.0 - eps);
			gTerm2 *= -gasRates[(int)Reaction.C1];

			double dmGas = (gTerm1 + gTerm2) * param.Axs * MW.CO2;  // A major term

			double term1b = pyro.gasYield * pyro.dryFeedIn.flow;    // A major term
			term1b *= pyro.dryGasOut.x[(int)Component.CO2] / Lbed;

			double[] mz = mDot[zIndex];
			double mzTar = mz[(int)Component.Tar];

			double term2b = pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);    // A major term
			term2b *= (mzTar * eps / U) * Arrhenius.wTarPyro[(int)Component.CO2];

			return term1b + term2b + dmGas;
		}


		// CH4 rate expression (dmCH4/dz)
		private double dmCH4(int zIndex, in double[] gasRates, double eps, double U)
		{
			double gTerm1 = -gasRates[(int)Reaction.SMR] * eps;

			double gTerm2 = (1.0 - eps) * (mChar / MW.Char);
			gTerm2 /= param.Axs * Lbed * (1.0 - eps);
			gTerm2 *= gasRates[(int)Reaction.C3];

			double dmGas = (gTerm1 + gTerm2) * param.Axs * MW.CH4;   // A major term

			double term1b = pyro.gasYield * pyro.dryFeedIn.flow;    // A major term
			term1b *= pyro.dryGasOut.x[(int)Component.CH4] / Lbed;

			double[] mz = mDot[zIndex];
			double mzTar = mz[(int)Component.Tar];

			double term2b = pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);    // A major term
			term2b *= (mzTar * eps / U) * Arrhenius.wTarPyro[(int)Component.CH4];

			return term1b + term2b + dmGas;
		}

		// H2 rate expression (dmH2/dz)
		private double dmH2(int zIndex, in double[] gasRates, double eps, double U)
		{
			double gTerm1 = eps * (3.0 * gasRates[(int)Reaction.SMR] + gasRates[(int)Reaction.WGS]);

			double gTerm2 = (1.0 - eps) * (mChar / MW.Char);
			gTerm2 /= param.Axs * Lbed * (1.0 - eps);
			gTerm2 *= gasRates[(int)Reaction.C2] - 2.0 * gasRates[(int)Reaction.C3];

			double dmGas = (gTerm1 + gTerm2) * param.Axs * MW.H2;  // A major term

			double term1b = pyro.gasYield * pyro.dryFeedIn.flow;    // A major term
			term1b *= pyro.dryGasOut.x[(int)Component.H2] / Lbed;

			double[] mz = mDot[zIndex];
			double mzTar = mz[(int)Component.Tar];

			double term2b = pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);    // A major term
			term2b *= (mzTar * eps / U) * Arrhenius.wTarPyro[(int)Component.H2];

			return term1b + term2b + dmGas;
		}

		// H2O rate expression (dmH2O/dz)
		private double dmH2O(int zIndex, in double[] gasRates, double eps)
		{
			double term1a = pyroGas.flowrate * pyroGas.x[(int)Component.H2O] / Lbed;

			double gTerm1 = -eps * (gasRates[(int)Reaction.SMR] + gasRates[(int)Reaction.WGS]);

			double gTerm2 = (1.0 - eps) * (mChar / MW.Char);
			gTerm2 /= param.Axs * Lbed * (1.0 - eps);
			gTerm2 *= -gasRates[(int)Reaction.C2];

			double term2a = (gTerm1 + gTerm2) * param.Axs * MW.H2O;

			return term1a + term2a;
		}


		// Tar rate expression (dmTar/dz)
		private double dmTar(int zIndex, double eps, double U)
		{
			double[] mz = mDot[zIndex];
			double term1 = pyro.tarYield * pyro.dryFeedIn.flow / Lbed;
			double term2 = -pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);
			term2 *= mz[(int)Component.Tar] * eps / U;
			return term1 + term2;
		}

	}
}

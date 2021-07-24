﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using BasicBFB.Model.Common;

namespace BasicBFB.Model
{
	public class ReactorBed
	{
		// --------------------------------------------------------------------------------
		//									PARAMETERS
		// --------------------------------------------------------------------------------
		public GasifierParams param;
		public Pyrolysis pyro;
		Stream pyroGas;

		List<double> zList = new List<double>(200);		
		List<double[]> mDotList = new List<double[]>(200);  // j = CO, CO2, CH4, H2, H2O, N2, O2, S, Tar

		public double[] mz0;		// double[Stream.numComp]

		public double mChar;
		public double Lbed;
		double epsAvg;

		const double TOL = 1.0e-5;      // Relative tolerance for mChar and Lbed convergence
		const double MAXITER = 50;


		// --------------------------------------------------------------------------------
		//									CONSTRUCTOR
		// --------------------------------------------------------------------------------
		public ReactorBed(GasifierParams param)
		{
			this.param = param;

			this.pyro = new Pyrolysis(param);
			pyro.pyrolize();

			// Initial guess for mChar holdup
			this.mChar = 0.05 * (pyro.dryFeedIn.flow * pyro.charYield);

			mixPyroGas();				// Instantiates and populates pyroGas
			setInitialConditions();		// Instantiates and populates mz0

			// Guess avgEps, then calculate Lbed using that
			epsAvg = 0.7;
			Lbed = param.L0 * (1.0 - param.epsMF) / (1.0 - epsAvg);
		}


		// --------------------------------------------------------------------------------
		//							MAIN CALCULATION ALGORITHM
		// --------------------------------------------------------------------------------
		public void solve()
		{
			zList.Clear();
			mDotList.Clear();
			double mCharErr = 1.0;
			double LbedErr = 1.0;
			bool isConverged = false;

			double[] zSpan = { 0.0, Lbed };
			OdeDel ode = dmDot;
			OdeSolver odeSolver = new OdeSolver();

			int iter0 = 0;
			int iter1 = 0;
			int iter2 = 0;

			while (!isConverged && (iter0 < MAXITER)) {
				iter0++;
				iter1 = 0;
				iter2 = 0;

				double mCharOld = mChar;
				while ((mCharErr > TOL) && (iter1 < MAXITER))
				{
					iter1++;
					(zList, mDotList) = odeSolver.rk45(ode, zSpan, mz0);
					updateMChar(zList, mDotList);
					mCharErr = Math.Abs(mChar - mCharOld) / mChar;
					mCharOld = mChar;
				}

				double LbedOld = Lbed;
				updateLbed(zList, mDotList);
				LbedErr = Math.Abs(Lbed - LbedOld) / Lbed;

				isConverged = ((mCharErr <= TOL) && (LbedErr <= TOL));

				while ((LbedErr > TOL) && (iter2 < MAXITER))
				{
					iter2++;
					LbedOld = Lbed;
					zSpan[1] = Lbed;
					(zList, mDotList) = odeSolver.rk45(ode, zSpan, mz0);
					updateLbed(zList, mDotList);
					LbedErr = Math.Abs(Lbed - LbedOld) / Lbed;
				}

				if (!isConverged)
				{
					updateMChar(zList, mDotList);
					mCharErr = Math.Abs(mChar - mCharOld) / mChar;
					isConverged = ( (mCharErr <= TOL) && (LbedErr <= TOL));
				}
			}

			if (isConverged)
			{
				Console.WriteLine("ReactorBed solve() converged after {0} cycles", iter0);
			} else
			{
				Console.WriteLine("ReactorBed solve() failed to converge!");
				Console.WriteLine($"mCharErr = {mCharErr:g4}, LbedErr = {LbedErr:g4}");
			}
		}



		// --------------------------------------------------------------------------------
		//							ODE FOR REACTION SYSTEM
		// --------------------------------------------------------------------------------
		private double[] dmDot(double z, in double[] mz)
		{
			double[] dmdz = new double[Stream.numComp];
			double[] gasRates = gasifyRates(z, mz);
			double eps = epsilon(z, mz);
			double U = gasU(mz);

			dmdz[(int)Component.CO] = dmCO(mz, gasRates, eps, U);
			dmdz[(int)Component.CO2] = dmCO2(mz, gasRates, eps, U);
			dmdz[(int)Component.CH4] = dmCH4(mz, gasRates, eps, U);
			dmdz[(int)Component.H2] = dmH2(mz, gasRates, eps, U);
			dmdz[(int)Component.H2O] = dmH2O(gasRates, eps);
			dmdz[(int)Component.Tar] = dmTar(mz, eps, U);
			dmdz[(int)Component.N2] = 0.0;
			dmdz[(int)Component.O2] = 0.0;
			dmdz[(int)Component.S] = 0.0;

			return dmdz;
		}


		// --------------------------------------------------------------------------------
		//							CHAR HOLDUP CALCULATIONS
		// --------------------------------------------------------------------------------

		// Calculates new value for mChar and updates this object's mChar parameter value
		private void updateMChar(in List<double> zs, in List<double[]> mdots)
		{
			double tauP = meanTauP(zs, mdots);
			Del3 rCsum = sumCharRates;
			double integral = Solver.integrate(rCsum, zs, mdots) / Lbed;
			double denom = (1.0 / tauP) + integral;
			mChar = pyro.charYield * pyro.dryFeedIn.flow / denom;
		}

		// d(tauP)/dz - derivative of particle residence time wrt z as a function of zIndex
		private double dtauP(double z, in double[] mz)
		{
			double eps = epsilon(z, mz);
			double U = gasU(mz);
			return (1.0 - eps) / U;
		}

		// Mean solids residence time tauP
		private double meanTauP(in List<double> zs, in List<double[]> mdots)
		{
			Del3 dtau = dtauP;
			return Solver.integrate(dtau, zs, mdots);
		}

		// Function used as integrand in mChar calculation method
		// Parameter z is not used, it's just so Solver.integrate can be reused
		private double sumCharRates(double z, in double[] mz)
		{
			double Tk = param.T + 273.15;       // T in Kelvin

			double pCO2 = partialP(Component.CO2, mz);
			double pH2O = partialP(Component.H2O, mz);
			double pH2 = partialP(Component.H2, mz);

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
			for (int i = 0; i < Stream.numComp; i++)
			{
				f[i] = pyro.dryGasOut.flowrate * pyro.dryGasOut.x[i]
						+ pyro.badGasOut.flowrate * pyro.badGasOut.x[i];
			}

			pyroGas.flowrate = pyro.dryGasOut.flowrate + pyro.badGasOut.flowrate;

			for (int i = 0; i < Stream.numComp; i++)
			{
				pyroGas.x[i] = f[i] / pyroGas.flowrate;
			}
		}

		// Add first point to z = 0, and add initial conditions to mDot
		private void setInitialConditions()
		{
			mz0 = new double[Stream.numComp];
			for (int i = 0; i < Stream.numComp; i++)
			{
				// Steam components are nonzero at z = 0, rest contribute nothing
				mz0[i] = param.steamIn.flowrate * param.steamIn.x[i];
			}
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
		//						BUBBLING BED PROPERTY CALCULATIONS
		// --------------------------------------------------------------------------------

		// Estimated bubble size db as a function of height (per zIndex)
		private double bubbleSize(double z, in double[] mz)
		{
			double db = 0.0;
			double deltaU = gasU(mz) - param.U0;
			double zTerm = z + 4.0 * Math.Sqrt((param.Axs / ((double)param.Nor)));
			db = 0.54 * Math.Pow(Const.g, -0.2) * Math.Pow(deltaU, 0.4);
			db *= Math.Pow(zTerm, 0.8);
			return db;
		}

		// Local bed voidage eps(z) based on zIndex
		private double epsilon(double z, in double[] mz)
		{
			double deltaU = gasU(mz) - param.U0;
			double dbTerm = 0.711 * Math.Sqrt(Const.g * bubbleSize(z, mz));
			double denom = 1.0 + deltaU / dbTerm;
			return 1.0 - (1.0 - param.epsMF) / denom;
		}

		// Average voidage of the whole fluidized bed
		private double avgEpsilon(in List<double> zs, in List<double[]> mDots)
		{
			Del3 eps = epsilon;
			double epsIntegral = Solver.integrate(eps, zs, mDots);
			return epsIntegral / Lbed;
		}

		// Calculate height of bubbling bed and set value of Lbed
		private void updateLbed(in List<double> zs, in List<double[]> mDots)
		{
			epsAvg = avgEpsilon(zs, mDots);
			Lbed = param.L0 * (1.0 - param.epsMF) / (1.0 - epsAvg);
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
		private double dmCO(in double[] mz, in double[] gasRates, double eps, double U)
		{
			double gTerm1 = eps * (gasRates[(int)Reaction.SMR] - gasRates[(int)Reaction.WGS]);
			
			double gTerm2 = (1.0 - eps) * (mChar / MW.Char);
			gTerm2 /= param.Axs * Lbed * (1.0 - eps);
			gTerm2 *= 2.0 * gasRates[(int)Reaction.C1] + gasRates[(int)Reaction.C2];
			
			double dmGas = (gTerm1 + gTerm2) * param.Axs * MW.CO;	// A major term

			double term1b = pyro.gasYield * pyro.dryFeedIn.flow;    // A major term
			term1b *= pyro.dryGasOut.x[(int)Component.CO] / Lbed;

			double mzTar = mz[(int)Component.Tar];

			double term2b = pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);    // A major term
			term2b *= (mzTar * eps / U) * Arrhenius.wTarPyro[(int)Component.CO];

			return term1b + term2b + dmGas;
		}


		// CO2 rate expression (dmCO2/dz)
		private double dmCO2(in double[] mz, in double[] gasRates, double eps, double U)
		{
			double gTerm1 = eps * gasRates[(int)Reaction.WGS];

			double gTerm2 = (1.0 - eps) * (mChar / MW.Char);
			gTerm2 /= param.Axs * Lbed * (1.0 - eps);
			gTerm2 *= -gasRates[(int)Reaction.C1];

			double dmGas = (gTerm1 + gTerm2) * param.Axs * MW.CO2;  // A major term

			double term1b = pyro.gasYield * pyro.dryFeedIn.flow;    // A major term
			term1b *= pyro.dryGasOut.x[(int)Component.CO2] / Lbed;

			double mzTar = mz[(int)Component.Tar];

			double term2b = pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);    // A major term
			term2b *= (mzTar * eps / U) * Arrhenius.wTarPyro[(int)Component.CO2];

			return term1b + term2b + dmGas;
		}


		// CH4 rate expression (dmCH4/dz)
		private double dmCH4(in double[] mz, in double[] gasRates, double eps, double U)
		{
			double gTerm1 = -gasRates[(int)Reaction.SMR] * eps;

			double gTerm2 = (1.0 - eps) * (mChar / MW.Char);
			gTerm2 /= param.Axs * Lbed * (1.0 - eps);
			gTerm2 *= gasRates[(int)Reaction.C3];

			double dmGas = (gTerm1 + gTerm2) * param.Axs * MW.CH4;   // A major term

			double term1b = pyro.gasYield * pyro.dryFeedIn.flow;    // A major term
			term1b *= pyro.dryGasOut.x[(int)Component.CH4] / Lbed;

			double mzTar = mz[(int)Component.Tar];

			double term2b = pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);    // A major term
			term2b *= (mzTar * eps / U) * Arrhenius.wTarPyro[(int)Component.CH4];

			return term1b + term2b + dmGas;
		}

		// H2 rate expression (dmH2/dz)
		private double dmH2(in double[] mz, in double[] gasRates, double eps, double U)
		{
			double gTerm1 = eps * (3.0 * gasRates[(int)Reaction.SMR] + gasRates[(int)Reaction.WGS]);

			double gTerm2 = (1.0 - eps) * (mChar / MW.Char);
			gTerm2 /= param.Axs * Lbed * (1.0 - eps);
			gTerm2 *= gasRates[(int)Reaction.C2] - 2.0 * gasRates[(int)Reaction.C3];

			double dmGas = (gTerm1 + gTerm2) * param.Axs * MW.H2;  // A major term

			double term1b = pyro.gasYield * pyro.dryFeedIn.flow;    // A major term
			term1b *= pyro.dryGasOut.x[(int)Component.H2] / Lbed;

			double mzTar = mz[(int)Component.Tar];

			double term2b = pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);    // A major term
			term2b *= (mzTar * eps / U) * Arrhenius.wTarPyro[(int)Component.H2];

			return term1b + term2b + dmGas;
		}

		// H2O rate expression (dmH2O/dz)
		private double dmH2O(in double[] gasRates, double eps)
		{
			double mH2O = param.feedIn.flow * param.feedIn.fracMoisture;
			double term1a = mH2O / Lbed;
			//double term1a = pyroGas.flowrate * pyroGas.x[(int)Component.H2O] / Lbed;

			double gTerm1 = -eps * (gasRates[(int)Reaction.SMR] + gasRates[(int)Reaction.WGS]);

			double gTerm2 = (1.0 - eps) * (mChar / MW.Char);
			gTerm2 /= param.Axs * Lbed * (1.0 - eps);
			gTerm2 *= -gasRates[(int)Reaction.C2];

			double term2a = (gTerm1 + gTerm2) * param.Axs * MW.H2O;

			return term1a + term2a;
		}


		// Tar rate expression (dmTar/dz)
		private double dmTar(in double[] mz, double eps, double U)
		{
			double term1 = pyro.tarYield * pyro.dryFeedIn.flow / Lbed;
			double term2 = -pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);
			term2 *= mz[(int)Component.Tar] * eps / U;
			return term1 + term2;
		}

	}
}
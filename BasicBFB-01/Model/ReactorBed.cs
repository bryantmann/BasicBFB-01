using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using BasicBFB.Model.Common;

namespace BasicBFB.Model
{
	public class ReactorBed : ReactorZone
	{
		// --------------------------------------------------------------------------------
		//									PARAMETERS
		// --------------------------------------------------------------------------------
		public Stream pyroGas { get; private set; }

		public double mChar { get; set; }				// Char holdup, kg
		public double Lbed { get; private set; }		// Height of expanded fluidized bed, m
		public double epsAvg { get; private set; }		// Average void fraction of bed
		public double mDotCharOut { get; private set; }	// Mass rate of char particles exiting the bed


		// --------------------------------------------------------------------------------
		//									CONSTRUCTOR
		// --------------------------------------------------------------------------------
		public ReactorBed(ref GasifierParams param)
		{
			this.param = param;

			this.pyro = new Pyrolysis(param);
			this.zList = new List<double>(200);
			this.mDotList = new List<double[]>(200);
			this.zSpan = new double[2];

			pyro.pyrolize();

			// Initial guess for mChar holdup
			this.mChar = 0.05 * (pyro.dryFeedIn.flow * pyro.charYield);
			this.mDotCharOut = 0.0;

			mixPyroGas();				// Instantiates and populates pyroGas
			setInitialConditions();		// Instantiates and populates mz0

			// Guess avgEps, then calculate Lbed using that
			epsAvg = 0.7;
			Lbed = param.L0 * (1.0 - param.epsMF) / (1.0 - epsAvg);
			this.effluent = new Effluent();
		}


		// --------------------------------------------------------------------------------
		//							MAIN CALCULATION ALGORITHMS
		// --------------------------------------------------------------------------------
		
		//public Effluent calcEffluent()
		//{
		//	solve();
		//	return effluent;
		//}

		public override void solve()
		{
			zList.Clear();
			mDotList.Clear();
			double mCharErr = 1.0;
			double LbedErr = 1.0;
			bool isConverged = false;

			zSpan[0] = 0.0;
			zSpan[1] = Lbed;
			OdeDel ode = dmDot;
			OdeSolver odeSolver = new OdeSolver();

			int iter0 = 0;
			int iter1 = 0;
			int iter2 = 0;

			while (!isConverged && (iter0 < Const.MAXITER)) {
				iter0++;
				iter1 = 0;
				iter2 = 0;

				double mCharOld = mChar;
				while ((mCharErr > Const.TOL) && (iter1 < Const.MAXITER))
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

				isConverged = ((mCharErr <= Const.TOL) && (LbedErr <= Const.TOL));

				while ((LbedErr > Const.TOL) && (iter2 < Const.MAXITER))
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
					isConverged = ( (mCharErr <= Const.TOL) && (LbedErr <= Const.TOL));
				}
			}

			if (isConverged)
			{
				Console.WriteLine("ReactorBed solve() converged after {0} cycles", iter0);
				int iMax = zList.Count - 1;     // Last index of zList and mDotList
				this.effluent = new Effluent(mDotList[iMax], zList[iMax], param, pyro, mDotCharOut);
			} else
			{
				Console.WriteLine("ReactorBed solve() failed to converge!");
				Console.WriteLine($"mCharErr = {mCharErr:g4}, LbedErr = {LbedErr:g4}");
				//this.effluent = new Effluent();
			}
		}



		// --------------------------------------------------------------------------------
		//							ODE for REACTION SYSTEM
		// --------------------------------------------------------------------------------
		
		protected override double[] dmDot(double z, in double[] mz)
		{
			double[] dmdz = new double[Stream.numComp];
			double[] gasRates = gasifyRates(z, mz);
			double eps = epsilon(z, mz);
			double U = gasU(mz);

			dmdz[(int)Component.CO] = dmCO(mz, gasRates, U, eps);
			dmdz[(int)Component.CO2] = dmCO2(mz, gasRates, U, eps);
			dmdz[(int)Component.CH4] = dmCH4(mz, gasRates, U, eps);
			dmdz[(int)Component.H2] = dmH2(mz, gasRates, U, eps);
			dmdz[(int)Component.H2O] = dmH2O(mz, gasRates, U, eps);
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
			//double tauP = meanTauP(zs, mdots);
			double tauP = param.tauP;
			Del3 rCsum = sumCharRates;
			double integral = Solver.integrate(rCsum, zs, mdots) / Lbed;
			double denom = (1.0 / tauP) + integral;
			mChar = pyro.charYield * pyro.dryFeedIn.flow / denom;

			mDotCharOut = mChar / tauP;
		}

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
		//						BUBBLING BED PROPERTY CALCULATIONS
		// --------------------------------------------------------------------------------

		// Estimated bubble size db as a function of height (per zIndex)
		private double bubbleSize(double z, in double[] mz)
		{
			double db = 0.0;
			double deltaU = gasU(mz) - param.Umf;
			double zTerm = z + 4.0 * Math.Sqrt((param.Axs / ((double)param.Nor)));
			db = 0.54 * Math.Pow(Const.g, -0.2) * Math.Pow(deltaU, 0.4);
			db *= Math.Pow(zTerm, 0.8);
			return db;
		}

		// Local bed voidage eps(z) based on zIndex
		internal override double epsilon(double z, in double[] mz)
		{
			double deltaU = gasU(mz) - param.Umf;
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
		internal override double dmCO2(in double[] mz, in double[] gasRates, double U, double eps)
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
		internal override double dmCH4(in double[] mz, in double[] gasRates, double U, double eps)
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
		internal override double dmH2(in double[] mz, in double[] gasRates, double U, double eps)
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
		internal override double dmH2O(in double[] mz, in double[] gasRates, double U, double eps)
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
		internal override double dmTar(in double[] mz, double U, double eps)
		{
			double term1 = pyro.tarYield * pyro.dryFeedIn.flow / Lbed;
			double term2 = -pyro.kPyros[3] * (1.0 - pyro.gammaTarInert);
			term2 *= mz[(int)Component.Tar] * eps / U;
			return term1 + term2;
		}

	}
}



//// d(tauP)/dz - derivative of particle residence time wrt z as a function of zIndex
//private double dtauP(double z, in double[] mz)
//{
//	double eps = epsilon(z, mz);
//	double U = gasU(mz);
//	return (1.0 - eps) / U;
//}

//// Mean solids residence time tauP
//private double meanTauP(in List<double> zs, in List<double[]> mdots)
//{
//	Del3 dtau = dtauP;
//	return Solver.integrate(dtau, zs, mdots);
//}

// Function used as integrand in mChar calculation method
// Parameter z is not used, it's just so Solver.integrate can be reused
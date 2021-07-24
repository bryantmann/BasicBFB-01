using System;
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

		public List<double> z = new List<double>();		
		public List<double[]> mDot = new List<double[]>();  // j = CO, CO2, CH4, H2, H2O, N2, O2, S, Tar

		public double mChar;
		public double Lbed;
		double epsAvg;


		// --------------------------------------------------------------------------------
		//									CONSTRUCTORS
		// --------------------------------------------------------------------------------
		public ReactorBed(GasifierParams param, double mCharGuess)
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
			double zTerm = z + 4.0 * Math.Sqrt((param.A / ((double)param.Nor)));
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



	}
}

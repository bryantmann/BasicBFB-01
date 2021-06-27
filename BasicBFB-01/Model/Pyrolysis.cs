using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicBFB
{
	public class Pyrolysis
	{
		public GasifierParams param { get; private set; }
		public double T { get; private set; }
		public Assay dryFeedIn { get; private set; }

		public double[] kPyros { get; private set; }    // All four [=] 1/sec
		public double[] phis { get; private set; }      // Dry organic tar mass fractions (n = 3)
		public double[] alphas { get; init; }           // Pyrolysis water mass fractions (n = 3)

		public double[] gasHrange { get; } = { 0.04, 0.105 };   // From gas stoichometry constraints
		public double rootTol = 1.0e-6;					// Absolute error for tar %H root 

		public double gasYield { get; private set; }	// Yields are wt/wt dry biomass
		public double tarYield { get; private set; }
		public double charYield { get; private set; }


		public Assay tarCHO { get; private set; }       // Assay on Wet basis
		public Assay dryGasCHO { get; private set;  }		// Stores gas CHO wt% values


		// Effluent streams
		public Stream dryGasOut { get; private set; }	// Dry, N2 free still
		public Stream badGasOut { get; private set; }		// Steam, N2, O2, H2S
		

		// TODO: Be sure to adjust dry feed rate down after cleanup and dryout



		// --------------------------------------------------------------------------
		// -------------------------- Constructors, etc -----------------------------
		// --------------------------------------------------------------------------

		public Pyrolysis(GasifierParams param)
		{
			this.param = param;			// Reference copy - single class instance for param
			this.T = param.T;
			this.dryFeedIn = param.feedIn.cleanedAndDried();
			this.dryFeedIn.name = param.feedIn.name + " (dry)";

			this.kPyros = Arrhenius.pyroRateConstants(param.T);
			this.gasYield = Arrhenius.pyroYield(kPyros, Phase.Gas);
			this.tarYield = Arrhenius.pyroYield(kPyros, Phase.Tar);
			this.charYield = Arrhenius.pyroYield(kPyros, Phase.Char);

			this.dryGasOut = new Stream(param.T, param.p);		// Mass fraction basis
			this.badGasOut = new Stream(param.T, param.p);
			
			this.tarCHO = new Assay("Tar");
			tarCHO.p = param.p;
			tarCHO.T = param.T;

			this.dryGasCHO = new Assay("Dry pyrolysis gas");
			dryGasCHO.p = param.p;
			dryGasCHO.T = param.T;


			phis = organicTarMassFrac();

			alphas = new double[3];				
			alphas[0] = 0.0;
			alphas[1] = 2 * MW.H / MW.H2O;
			alphas[2] = MW.O / MW.H2O;
		}


		// --------------------------------------------------------------------------
		// -------------------------------- Methods ---------------------------------
		// --------------------------------------------------------------------------

		// True if no converged ok, false if exception was thrown
		public bool pyrolize()
		{
			// Step 1: Calculate tar %H using root finder algorithm.
			Del1 handler = tarBalance;
			double gasH = 0.0;
			try
			{
				gasH = Solver.ridder(tarBalance, gasHrange, rootTol);
			}
			catch (NotConvergedException e)
			{
				Console.WriteLine(e.Message);
				return false;
			}
			catch (ArgumentException e)
			{
				Console.WriteLine(e.Message);
				return false;
			}

			// Step 2: Calculate and store the tar wt% CHO in the tarOut assay
			tarCHO.w[1] = tarFracFromGas(Element.H, gasH);
			tarCHO.w[0] = tarFracFromTar(Element.C, gasH);
			tarCHO.w[2] = tarFracFromTar(Element.O, gasH);

			// Step 3: using tar composition, populate this.dryGasCHO.w results
			dryGasCHO.w[1] = gasH;
			dryGasCHO.w[0] = gasFracFromTar(Element.C, tarCHO.w[0]);
			dryGasCHO.w[2] = gasFracFromTar(Element.O, tarCHO.w[2]);

			// Step 4: Populate dryGasOut stream components from aboce using MATH!
			dryGasOut.setX(gasComposition());
			if (!dryGasOut.isNormalized())
			{
				string objName = "dryGasOut";
				double sum = dryGasOut.sumX();
				string message = $"Stream {objName} has sum(x[i]) = {sum}, above the normalization tolerance";
				Console.WriteLine(message);
				throw new UnexpectedValueException(objName, sum, message);

			}
			
			// Step 5: Populate badGasOut

			// Step 6: Adjust flowrates, also with MATH!

			// Step 7: Profit!

			return true;
		}


		// Given a value for wHgas, calculate tar CHO values.  Target is for sum = 1.000
		// This method returns the difference of that sum from 1.0 (target zero difference)
		// Used by a delegate in ridder root finding method
		private double tarBalance(double wHgas)
		{
			double[] wTar = new double[3];

			// Calculate tar hydrogen from gas H
			wTar[1] = tarFracFromGas(Element.H, wHgas);

			// Next calculate tar %C and %O from tar %H
			wTar[0] = tarFracFromTar(Element.C, wTar[1]);
			wTar[2] = tarFracFromTar(Element.O, wTar[1]);

			// Add them up and subtract from 1.0 to return result of balance
			double sum = 0.0;
			foreach (double w in wTar)
			{
				sum += w;
			}

			return 1.0 - sum;
		}



		// TODO: Pyrolysis water method

		private double[] organicTarMassFrac()
		{
			double[] xDryTar = new double[3];
			xDryTar[0] = 1.14 * dryFeedIn.w[0];
			xDryTar[1] = 1.13 * dryFeedIn.w[1];
			xDryTar[2] = 1.0 - xDryTar[0] - xDryTar[1];
			return xDryTar;
		}


		// Calculates the tar mass fraction of the specified element
		// based on its mass fraction in the tar (passed as argument wG)
		private double tarFracFromGas(Element element, double wG)
		{
			double wB = dryFeedIn.w[(int)element];
			double wC = (int)element == 0 ? 1.0 : 0.0;

			double wT = (kPyros[0] + kPyros[1] + kPyros[2]) * wB;
			wT = wT - kPyros[0] * wG - kPyros[2] * wC;
			wT /= kPyros[1];

			return wT;
		}


		private double tarFracFromTar(Element element, double wHtar)
		{
			// Element specific values are singular nouns
			double phi = phis[(int)element];
			double alpha = alphas[(int)element];

			double wT = phi - (alpha - phi) * (wHtar - phis[1]) / (alphas[1] - phis[1]);
			return wT;
		}


		private double gasFracFromTar(Element element, double wT)
		{
			double S = kPyros[0] + kPyros[1] + kPyros[2];
			double wB = dryFeedIn.w[(int)element];
			double wC = (int)element == 0 ? 1.0 : 0.0;

			double wG = S * wB - kPyros[1] * wT - kPyros[2] * wC;
			wG /= kPyros[0];

			return wG;
		}


		// Calculates CO, CO2, CH4 and H2 in dry, nitrogen-free pyrolysis gas
		// Uses same indexing as Stream and the Component enumeration
		private double[] gasComposition()
		{
			double[] y = new double[10];        // Will only fill the first 4 spots

			// Experiments indicate H2 from pyrolysis is so small it can be neglected'
			y[3] = 0.0;							// y[3] is H2
			y[2] = 4.0 * dryGasCHO.w[1];           // w[1] is wt frac H

			// Next we can calc CO2.  I'll do it piecewise as usual
			double alpha1 = (12.0 / 16.0) * y[3];
			double beta1 = (924.0 / 28.0) * dryGasCHO.w[2];		// w[2] is O
			double gamma1 = (3.0 / 11.0 - 6.0 / 7.0);	
			y[1] = (dryGasCHO.w[0] - alpha1 - beta1) / gamma1;     // y[2] is CO2, w[0] is C

			// CO is last
			double a = (12.0 / 16.0) * y[3];
			double b = (12.0 / 44.0) * y[1];
			double c = 12.0 / 28.0;
			y[0] = (dryGasCHO.w[0] - a - b) / c;

			return y;
		}

	

	}
}

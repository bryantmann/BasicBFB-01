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

		public double[] gasHrange { get; } = { 0.04, 0.105 };	// From gas stoichometry constraints
		
		public double gasYield { get; private set; }	// Yields are wt/wt dry biomass
		public double tarYield { get; private set; }
		public double charYield { get; private set; }



		public Stream dryGasOut { get; private set; }	
		public Stream badGasOut { get; private set; }   // Steam, N2, O2, H2S
		public Assay tarOut { get; private set; }		// Assay on Wet basis

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
			
			this.tarOut = new Assay();
			tarOut.p = param.p;
			tarOut.T = param.T;

			phis = organicTarMassFrac();

			alphas = new double[3];				
			alphas[0] = 0.0;
			alphas[1] = 2 * MW.H / MW.H2O;
			alphas[2] = MW.O / MW.H2O;
		}


		// --------------------------------------------------------------------------
		// -------------------------------- Methods ---------------------------------
		// --------------------------------------------------------------------------


		public double[] organicTarMassFrac()
		{
			double[] xDryTar = new double[3];
			xDryTar[0] = 1.14 * dryFeedIn.w[0];
			xDryTar[1] = 1.13 * dryFeedIn.w[1];
			xDryTar[2] = 1.0 - xDryTar[0] - xDryTar[1];
			return xDryTar;
		}


		// Calculates the tar mass fraction of the specified element
		// based on its mass fraction in the tar (passed as argument wG)
		double tarFracFromGas(Element element, double wG)
		{
			double wB = dryFeedIn.w[(int)element];
			double wC = (int)element == 0 ? 1.0 : 0.0;

			double wT = (kPyros[0] + kPyros[1] + kPyros[2]) * wB;
			wT = wT - kPyros[0] * wG - kPyros[2] * wC;
			wT /= kPyros[1];

			return wT;
		}


		double tarFracFronTar(Element element, double wHtar)
		{
			// Element specific values are singular nouns
			double phi = phis[(int)element];
			double alpha = alphas[(int)element];

			double wT = phi - (alpha - phi) * (wHtar - phis[1]) / (alphas[1] - phis[1]);
			return wT;
		}


		// MARK: - The Main Event ---------------------------------------------------------

		// For root finding algorithm use
		// Given a value for wHgas, calculate tar CHO values.  Target is for sum = 1.000
		// This method returns the difference of that sum from 1.0 (target zero difference)
		double tarBalance(double wHgas)
		{
			double[] wTar = new double[3];

			// Calculate tar hydrogen from gas H
			wTar[1] = tarFracFromGas(Element.H, wHgas);

			// Next calculate tar %C and %O from tar %H
			wTar[0] = tarFracFronTar(Element.C, wTar[1]);
			wTar[2] = tarFracFronTar(Element.O, wTar[1]);

			// Add them up and subtract from 1.0 to return result of balance
			double sum = 0.0;
			foreach (double w in wTar)
			{
				sum += w;
			}

			return 1.0 - sum;
		}






		// --------------------------------------------------------------------------------
		// TODO - Get rid of pyroGasH(xC) and pyroGasC(xH), they suck.

		// Returns mass fraction of hydrogen in cracked gas for a given value of xC
		// From mass elemental balance.  
		public double pyroGasH(double xC)
		{
			double[] phi = organicTarMassFrac();
			double S = kPyros[0] + kPyros[1] + kPyros[2];
			double[] alpha = new double[3];			// Mass fractions of C, H, O in H2O
			alpha[0] = 0.0;
			alpha[2] = MW.O / MW.H2O;
			alpha[1] = 1.0 - alpha[2];

			double[] beta = new double[3];          // Mass fractions in dry, clean feed
			Array.Copy(dryFeedIn.w, 0, beta, 0, 3);

			double eta = (alpha[2] - phi[2]) / phi[0];

			// Temporary storage array to simplify code
			double[] terms = new double[6];
			terms[0] = 1.0;
			terms[1] = -xC;
			terms[2] = -beta[2] * S / kPyros[0];
			terms[3] = phi[2] * kPyros[1] / kPyros[0];

			terms[4] = beta[0] * S - kPyros[0] * xC - kPyros[2];
			terms[4] *= -eta / kPyros[1];

			terms[5] = phi[0] * kPyros[1] / kPyros[0];
				
			double xH = 0.0;
			foreach (double t in terms)
			{
				xH += t;
			}

			return xH;
		}


		// Returns mass fraction of carbon in cracked gas for a given value of xH
		// From mass elemental balance.  Will cross-check these two methods
		public double pyroGasC(double xH)
		{
			double[] phi = organicTarMassFrac();
			double S = kPyros[0] + kPyros[1] + kPyros[2];
			double[] alpha = new double[3];         // Mass fractions of C, H, O in H2O
			alpha[0] = 0.0;
			alpha[2] = MW.O / MW.H2O;
			alpha[1] = 1.0 - alpha[2];

			double[] beta = new double[3];          // Mass fractions in dry, clean feed
			//Array.Copy(dryFeedIn.w, 0, beta, 0, 3);

			double gamma = (alpha[2] - phi[2]) / (alpha[1] - phi[1]);

			double[] terms = new double[6];
			terms[0] = 1.0;
			terms[1] = -xH;
			terms[2] = -phi[2] * kPyros[1] / kPyros[0];
			terms[3] = -beta[2] * S / kPyros[0];

			terms[4] = beta[1] * S - kPyros[0] * xH;
			terms[4] *= gamma / kPyros[0];

			terms[5] = -gamma * phi[1] * kPyros[1] / kPyros[0];

			double xC = 0.0;
			foreach (double t in terms)
			{
				xC += t;
			}

			return xC;
		}

	}
}

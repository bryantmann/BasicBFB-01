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

		public double[] kPyros { get; private set; }	// All four [=] 1/sec
		public double gasYield { get; private set; }	// Yields are wt/wt dry biomass
		public double tarYield { get; private set; }
		public double charYield { get; private set; }


		public Stream dryGasOut { get; private set; }	
		public Stream wetGasOut { get; private set; }   // Steam, N2, O2, H2S
		public Assay tarOut { get; private set; }		// Assay on Wet basis
		public Assay charOut { get; private set; }		// Char assay without ash
		public Assay ashesOut { get; private set; }

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
			this.wetGasOut = new Stream(param.T, param.p);
			
			this.tarOut = new Assay();
			tarOut.p = param.p;
			tarOut.T = param.T;

			this.charOut = tarOut.Clone();
			this.ashesOut = tarOut.Clone();
		}


		// --------------------------------------------------------------------------
		// -------------------------------- Methods ---------------------------------
		// --------------------------------------------------------------------------

		// Dry, organic Tar composition
		public Assay gasFromOrganicTar()
		{
			double[] xDryTar = new double[3];
			xDryTar[0] = 1.14 * dryFeedIn.w[0];
			xDryTar[1] = 1.13 * dryFeedIn.w[1];
			xDryTar[2] = 1.0 - xDryTar[0] - xDryTar[1];
			
			Assay wfDryTar = new Assay(xDryTar);
			wfDryTar.T = this.T;
			wfDryTar.p = this.T;
			
			return wfDryTar;
		}

		public double[] organicTarMassFrac()
		{
			double[] xDryTar = new double[3];
			xDryTar[0] = 1.14 * dryFeedIn.w[0];
			xDryTar[1] = 1.13 * dryFeedIn.w[1];
			xDryTar[2] = 1.0 - xDryTar[0] - xDryTar[1];
			return xDryTar;
		}


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

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


		//// For a given potential solution for wt%H in wet tar, determine %C, %O
		//public Assay pyroTarAssay()
		//{
		//	Assay wfWetTar = new Assay();
		//	wfWetTar.T = this.T;
		//	wfWetTar.p = this.p;
		//	wfWetTar.flow = dryFeedIn.flow * tarYield;
		//	// More code goes here

		//	return wfWetTar;
		//}


		// Returns mass fraction of hydrogen in cracked gas for a given value of xC
		// From mass elemental balance.  
		public double pyroGasH(double xC)
		{
			Assay dryTarGas = gasFromOrganicTar();
			double sumK = kPyros[0] + kPyros[1] + kPyros[2];
			double xH = 0.0;

			// Break giant formula up into pieces
			double[] termsA = new double[5];
			double[] termsB = new double[3];
			double[] factors = new double[3];

			termsA[0] = 1.0;
			termsA[1] = -xC;
			termsA[2] = dryTarGas.w[2] * kPyros[1] / kPyros[0];
			termsA[3] = -sumK * dryFeedIn.w[2] / kPyros[0];

			// The last term is large and can be grouped into factors and subterms
			termsB[0] = sumK * dryFeedIn.w[0];
			termsB[1] = -kPyros[0] * xC;
			termsB[2] = kPyros[2];

			double sumBTerms = 0.0;
			foreach (double tB in termsB)
			{
				sumBTerms += tB;
			}

			factors[0] = kPyros[1] / kPyros[1];
			factors[1] = (MW.O / MW.H2O) - dryTarGas.w[2];
			factors[2] = 1.0 - (1.0 / (kPyros[1] * dryTarGas.w[1])) * sumBTerms;

			// Now combine everything
			termsA[4]= factors[0] * factors[1] * factors[2];

			foreach (double tA in termsA)
			{
				xH += tA;
			}
			return xH;
		}


		// Returns mass fraction of carbon in cracked gas for a given value of xH
		// From mass elemental balance.  Will cross-check these two methods
		public double pyroGasC(double xH)
		{
			Assay dryTarGas = gasFromOrganicTar();
			double sumK = kPyros[0] + kPyros[1] + kPyros[2];
			double xC = 0.0;

			// Break this other giant formula into different pieces
			double[] termsA = new double[5];
			double[] termsB = new double[3];
			double[] factors = new double[2];

			termsA[0] = 1.0;
			termsA[1] = -xH;
			termsA[2] = dryTarGas.w[2] * kPyros[1] / kPyros[0];
			termsA[3] = -sumK * dryFeedIn.w[2] / kPyros[0];

			// Last term is huge, break into pieces
			termsB[0] = sumK * dryFeedIn.w[0];
			termsB[1] = -xH * kPyros[0];
			termsB[2] = -kPyros[1] * dryTarGas.w[0];

			factors[0] = 0.0;
			foreach (double tB in termsB)
			{
				factors[0] += tB;
			}

			factors[1] = (MW.O / MW.H2O) - dryTarGas.w[2];
			factors[1] /= kPyros[0] * ((MW.H / MW.H2O) - dryTarGas.w[0]);
			termsA[4] = factors[0] + factors[1];

			foreach (double tA in termsA)
			{
				xC += tA;
			}

			return xC;
		}

	}
}

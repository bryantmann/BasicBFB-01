using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicBFB
{
	// Arrhenius expression used to calculate rate constants vs temperature
	// This version of the model will mainly use this for pyrolysis reactions
	static public class Arrhenius
	{
		// Calculates the rate constants k1-k4 characterizing biomass pyrolysis
		// Results are returned in an array of length 4 (with k1 at index 0)

		const int numRxns = 4;

		static public double[] pyroRateConstants(double T)
		{
			double[] kVals = new double[numRxns];
			double kelvinT = T + 273.15;              // Requires T provided in °C
			
			for (int i = 0; i < numRxns; i++)
			{
				double exponent = PyroParams.Ea[i] * 1000.0;	// convert Ea [=] kJ/mol to J/mol
				exponent /= (Const.Rgas * kelvinT);
				kVals[i] = PyroParams.A0[i] * Math.Exp(-exponent);
			}

			return kVals;
		}


		// The following method calculates yields from the first set of pyro reactions
		// The Phase enumeration type is provided as one of the arguments
		static public double pyroYield(double T, Phase phi)
		{
			double yield = 0.0;
			double[] kPyro = pyroRateConstants(T);
			double kx = kPyro[(int)phi];	// Value depends on choice of Phase value
											// Enum raw values already chosen to line up

			double ksum = kPyro[0] + kPyro[1] + kPyro[2];
			yield = kx / ksum;

			return yield;
		}


		// Same thing but rate constants are passed as an argument
		static public double pyroYield(double[] kPyro, Phase phi)
		{
			double yield = 0.0;
			double kx = kPyro[(int)phi];    // Value depends on choice of Phase value
											// Enum raw values already chosen to line up

			double ksum = kPyro[0] + kPyro[1] + kPyro[2];
			yield = kx / ksum;

			return yield;
		}


		// Method to estimate the dry elemental composition of
		// "organic tar" that has emulsified water mixed in based on feed assay
		// Note the input values correspond to FULLY DRY biomass (moisture subtracted)
		static public Assay wfDryTar(Assay wfBiomass)
		{
			// Start by making a copy 
			Assay wfAdjusted = wfBiomass.Clone();

			

			
			// Simplify: Since N, S and metals are present in small amounts vs. CHO
			// So zero out all other element wt fractions and renormalize
			for (int i = 3; i < Assay.numElements; i++)
			{
				wfAdjusted.w[i] = 0;
			}

			

			return wfAdjusted;
		}


	}


	// The readonly struct PyroParams is simply a place to store constants
	// Indices for k array correspond to:
	//		0 - Gas formation pyrolysis reaction (k1 in published paper)
	//		1 - Tar formation pyrolysis reaction 
	//		2 - Char formation 
	//		3 - Secondary gas formation reaction from tar pyrolysis
	public readonly struct PyroParams
	{
		// Arrhenius parameters for two-stage pyrolysis mechanism
		public static double[] A0 = { 1.30e8, 2.00e8, 1.08e7, 1.0e5 }; // 1/s
		public static double[] Ea = { 140.0, 133.0, 121.0, 93.3 };        // kJ/mol
		public static double[] dHrxns = { 64.0, 64.0, 64.0, -42.0 };       // kJ/kg

		// Model assumes a fixed gas composition for gas from tar pyrolysis (!)
		// Temperature and pressure values aren't being captured in this Stream object though
		// Units of x are mass fraction
		public static Stream tarPyroGas = new Stream(new double[] { 0.7222, 0.1422, 0.1133, 0.0222 }); 
	}
}

/*
 *		Component fraction order in Stream objects' x array:
 *		  INDEX		  SPECIES
 *			0			CO
 *			1			CO2
 *			2			CH4		(all hydrocarbons lumped in with CH4)
 *			3			H2
 *			4			H2O
 *			5			N2		(may not be used much)
 *			6			O2		(may not be used much)
 *			7			H2S		(may not be used much)
 *			8			Tar		(may not be used much)
 *			9			Char	(may not be used much)
 */
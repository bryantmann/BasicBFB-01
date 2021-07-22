using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using BasicBFB.Model;


namespace BasicBFB.Testing
{
	// This static class generates input data used to run gasifier model tests
	// Source: https://doi.org/10.1021/acs.energyfuels.6b03161
	public static class ExampleData
	{
		
		// Creates an Assay object for "sawdust" biomass feed
		public static Assay sawdustAssay()
		{
			double[] w = new double[5] { 0.425, 0.063, 0.51, 0.002, 0.0 };
			Assay assay = new Assay(w);
			assay.fracMoisture = 0.085;
			assay.fracAsh = 0.012;
			assay.fracVolatiles = 0.774;
			assay.fracFixedCarbon = 0.129;
			assay.feedParticleSize = 500.0;     // microns
			return assay;
		}


		// 100% H2O at 100 C and 1 atm (T, P arbitrary chosen, can change)
		// Note the flowrate here is 0.0
		public static Stream pureSteam()
		{
			double[] w = new double[8];
			for (int i = 0; i < w.Length; i++)
			{
				w[i] = 0.0;
			}

			w[4] = 1.0;     // Index for H2O
			Stream steam = new Stream(w, temp: 100.0, press: 1.01325);
			return steam;
		}

		
		public static GasifierParams gasifierParams()
		{
			GasifierParams p = new GasifierParams();

			Assay sawdust = sawdustAssay();
			sawdust.name = "Sawdust";

			Stream steam = pureSteam();

			p.setVesselParams(T: 750.0, p: 1.115, dBed: 0.15, dFree: 0.15, Ltot: 1.2);
			
			p.setBedParams(U0: 0.25, L0: 0.32, mSand: 8.0, dSand: 250.0, rhoSand: 2600.0);

			p.setFeedParams(source: "Sawdust", feedRate: 3.75,
							dFeed: 500.0, feedIn: sawdust, stmIn: steam, sbRatio: 2.0);

			p.setupStreams();
			return p;
		}


	}
}

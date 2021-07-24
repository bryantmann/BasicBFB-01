using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using BasicBFB.Model;
using BasicBFB.Model.Common;

namespace BasicBFB.Testing
{
	// This static class generates input data used to run gasifier model tests
	// Source: https://doi.org/10.1021/acs.energyfuels.6b03161
	public static class ExampleData
	{
		private static double p0 = 1.115;   // bara
		private static double T0 = 750.0;	// degC


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
			double[] y = new double[Stream.numComp];
			for (int i = 0; i < w.Length; i++)
			{
				y[i] = 0.0;
			}

			y[(int)Component.H2O] = 0.9;     // Index for H2O
			y[(int)Component.N2] = 0.1;

			Stream steam = new Stream(w, temp: T0, press: p0, true);
			return steam;
		}

		
		public static GasifierParams gasifierParams()
		{
			GasifierParams p = new GasifierParams();

			Assay sawdust = sawdustAssay();
			sawdust.name = "Sawdust";

			p.setVesselParams(T: T0, p: p0, dBed: 0.15, dFree: 0.15, Ltot: 1.2);
			
			p.setBedParams(U0: 0.25, L0: 0.32, Nor: 50, mSand: 8.0, dSand: 250.0, rhoSand: 2600.0);

			p.setFeedParams(source: "Sawdust", feedRate: 3.75,
							dFeed: 500.0, feedIn: sawdust, stmIn: steam, sbRatio: 2.0);

			Stream steam = pureSteam();

			p.setupStreams();
			return p;
		}


	}
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using BasicBFB.Model.Common;

namespace BasicBFB.Model
{
	// Stores and communicates global conditions characterizing the system
	public class GasifierParams
	{
		// --------------------------------------------------------------------------
		// Vessel parameters and general conditions
		public double T { get; private set; }				// Reactor temp, °C (assumed const)
		public double p { get; private set; }				// Reactor pressure, bara (assumed const)
		public double dBed { get; private set; }            // Bed diameter in meters
		public double dFree { get; private set; }			// Freeboard diameter, m
		public double Ltot { get; private set; }            // Total gasifier vessel height, m
		public double Axs									// Reactor xsec area, units: m2
		{
			get
			{
				return Const.pi * dBed * dBed / 4.0;
			}
		}

		// --------------------------------------------------------------------------
		// Fluidized bed parameters
		public double L0 { get; private set; }              // Bed height at min U0 gas velocity, m
		public double U0 { get; private set; } = 0.25;		// Superficial gas velocity at inlet, m/s
		public int Nor { get; private set; }				// Number of orifices in distributor plate
		public double mSand { get; private set; }			// Inventory (kg) of sand in reactor
		public double dSand { get; private set; }			// Avg. sand particle diameter, µm
		public double rhoSand { get; private set; }         // Sand particle density, kg/m3
		public double epsMF { get; } = 0.41408;				
		public double Umf { get; private set; } = 0.025;	
		public double tauP { get; private set; }   // Mean solids restime (s)
					// Calculated from alpha using Umf, U0, mB and sbr

		// Alpha is optimized so mDot(char) = 0 at z = Ltot (freeboard effluent)
		public double alpha { get; set; } = 0.9;            // Char conversion in bed


		// --------------------------------------------------------------------------
		// Incoming biomass feed and gasifying agent		
		public string bioSource { get; private set; }		// Brief descriptor of feed type
		public double bioFeedRate { get; private set; }		// Units of kg/s
		public double dFeed { get; private set; }			// Avg. biomass particle diameter, µm
		public Assay feedIn { get; private set; }        // Feed assay properties
		public Stream steamIn { get; private set; }			// Gasifying stream properties
		public double sbr { get; private set; }             // Steam to biomass ratio, wt/wt
		public double[] ySteam { get; private set; }

		// --------------------------------------------------------------------------


		// --------------------------------------------------------------------------
		//									METHODS
		// --------------------------------------------------------------------------
		public void setVesselParams(double T, double p, double dBed, 
									double dFree, double Ltot)
		{
			this.T = T;
			this.p = p;
			this.dBed = dBed;
			this.dFree = dFree;
			this.Ltot = Ltot;
		}


		public void setBedParams(double U0, double L0, int Nor,
								 double mSand, double dSand, double rhoSand)
		{ 
			this.U0 = U0;
			this.L0 = L0;
			this.Nor = Nor;
			this.mSand = mSand;
			this.dSand = dSand;
			this.rhoSand = rhoSand;
		}


		public void setFeedParams(string source, double feedRate, double dFeed,
								  Assay feedIn, double[] ySteam, double sbRatio)
		{
			this.bioSource = source;
			this.bioFeedRate = feedRate;
			this.dFeed = dFeed;
			this.feedIn = feedIn;
			this.sbr = sbRatio;
			this.ySteam = ySteam;
		}


		public void setupStreams()
		{
			feedIn.flow = bioFeedRate;
			feedIn.p = this.p;
			feedIn.T = this.T;

			double[] nDotStm = new double[Stream.numComp];
			int jH2O = (int)Component.H2O;
			nDotStm[jH2O] = this.sbr * bioFeedRate / MW.H2O;
			for (int j = 0; j < Stream.numComp; j++)
			{
				nDotStm[j] = nDotStm[jH2O] * ySteam[j] / ySteam[jH2O];
			}
			steamIn = new Stream(T, p, false);
			steamIn.nDot = nDotStm;
		}


		public void calcFluidBedParams(bool calcU0 = true)
		{
			// Local units are cm, g, min or s - converted back when setting property
			double rhoGas = steamIn.rho;		// [=] g/cm3
			double rhoP = rhoSand * 1.0e-3;     // [=] g/cm3
			double dP = dSand * 1.0e-4;         // [=] cm
			double muGas = steamVis(T);			// [=] g/cm.s
			double g = Const.g * 100.0;			// cm/s2

			// Calculate Umf
			double Ar = rhoGas * dP * dP * dP * (rhoP - rhoGas) * g;
			Ar /= (muGas * muGas);
			// Quadratic formula - solution is positive root
			double aQ = 24.5;
			double bQ = 1650.0;
			double cQ = -Ar;
			double term2 = Math.Sqrt(bQ * bQ - 4.0 * aQ * cQ);
			double Re_mf = ((-bQ + term2) / (2.0 * aQ));

			double U_mf = (Re_mf * muGas) / (rhoGas * dP);  // cm/sec
			this.Umf = U_mf / 100.0;		// m/s

			if (calcU0)
			{
				double mDotGas = steamIn.mDotTotal * 1000.0;    // g/sec
				double Acm2 = Axs * 1.0e4;
				this.U0 = mDotGas / (rhoGas * Acm2);
			}

			// Calculate tauP (C.E. Agu et al, Fuel 253(2019) 1414-1423)
			double mDot_B = bioFeedRate * 60.0;     // kg/min
			double mP = mSand;                      // kg
			double Uratio = U0 / Umf;
			double xb = 4055.0 * (mDot_B / mP) * Math.Pow((Uratio), -0.185);
			xb = Math.Pow(xb, 1.385);

			// Both td and te in minutes
			double td = 11.35 * Math.Pow(xb, 0.028) * Math.Pow(Uratio, -0.3);
			double te = 67.58 * Math.Pow(xb, 0.278) * Math.Pow(Uratio, -0.185);

			tauP = (td + alpha * (te - td)) * 60.0;
		}

		// Viscosity of pure steam at p = 1.115 bara (500-800C data)
		// Units: Poise [=] g/cm.s		T [=] °C
		public double steamVis(double T)
		{
			const double slope = 3.935e-7;
			const double yInt = 8.968e-5;
			double mu = slope * T + yInt;
			return mu;
		}


		//public void setFeedParams(string source, double feedRate, double dFeed,
		//						  Assay feedIn, Stream stmIn, double sbRatio)
		//{
		//	this.bioSource = source;
		//	this.bioFeedRate = feedRate;
		//	this.dFeed = dFeed;
		//	this.feedIn = feedIn;
		//	this.steamIn = stmIn;
		//	this.sbr = sbRatio;
		//}

		// Used to set flowrates, p and T 
		//public void setupStreams()
		//{
		//	feedIn.flow = bioFeedRate;
		//	feedIn.p = this.p;
		//	feedIn.T = this.T;
		//	steamIn.p = this.p;		// For this model everything's constant T, P
		//	steamIn.T = this.T;		// I'm keeping these here for a more complex model

		//	if (steamIn.x[(int)Component.H2O] <= 0)
		//	{
		//		return;		// Don't set rate for steamIn if there's zero wt% water
		//	}

		//	if (!steamIn.isMolar)
		//	{
		//		steamIn.flowrate = sbr * bioFeedRate / steamIn.x[4];        // H2O is index 4 for Strean
		//	} else
		//	{
		//		// When steamIn.isMolar == true, both x and flowrate are molar
		//		// Calculate mass rate of steam from sbr, then convert to moles
		//		double wH2O = steamIn.x[4] * MW.H2O / steamIn.avgMW;
		//		double mDotTot = sbr * bioFeedRate / wH2O;
		//		steamIn.flowrate = mDotTot;		// Mass flow rate
		//		//steamIn.flowrate = mDotTot / steamIn.avgMW;
		//	}
		//}

	}
}

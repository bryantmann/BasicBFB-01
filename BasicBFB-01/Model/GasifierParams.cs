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

		// Vessel parameters and general conditions
		public double T { get; private set; }				// Reactor temp, °C (assumed const)
		public double p { get; private set; }				// Reactor pressure, bara (assumed const)
		public double dBed { get; private set; }            // Bed diameter in meters
		public double dFree { get; private set; }			// Freeboard diameter, m
		public double Ltot { get; private set; }            // Total gasifier vessel height, m
		
		
		// Fluidized bed parameters
		public double L0 { get; private set; }				// Bed height at min U0 gas velocity, m
		public double U0 { get; private set; }				// Minimum gas velocity, m/s 
		public int Nor { get; private set; }				// Number of orifices in distributor plate
		public double mSand { get; private set; }			// Inventory (kg) of sand in reactor
		public double dSand { get; private set; }			// Avg. sand particle diameter, µm
		public double rhoSand { get; private set; }         // Sand particle density, kg/m3
		public double epsMF = 0.41408;


		// Incoming biomass feed and gasifying agent		
		public string bioSource { get; private set; }		// Brief descriptor of feed type
		public double bioFeedRate { get; private set; }		// Units of kg/s
		public double dFeed { get; private set; }			// Avg. biomass particle diameter, µm
		public Assay feedIn { get; private set; }        // Feed assay properties
		public Stream steamIn { get; private set; }			// Gasifying stream properties
		public double sbr { get; private set; }             // Steam to biomass ratio, wt/wt

		public double Axs
		{
			get
			{
				return Const.pi * dBed * dBed / 4.0;
			}
		}


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
								  Assay feedIn, Stream stmIn, double sbRatio)
		{
			this.bioSource = source;
			this.bioFeedRate = feedRate;
			this.dFeed = dFeed;
			this.feedIn = feedIn;
			this.steamIn = stmIn;
			this.sbr = sbRatio;
		}


		// Used to set flowrates, p and T 
		public void setupStreams()
		{
			feedIn.flow = bioFeedRate;
			feedIn.p = this.p;
			feedIn.T = this.T;
			steamIn.p = this.p;		// For this model everything's constant T, P
			steamIn.T = this.T;		// I'm keeping these here for a more complex model
			
			if (steamIn.x[(int)Component.H2O] <= 0)
			{
				return;		// Don't set rate for steamIn if there's zero wt% water
			}

			if (!steamIn.isMolar)
			{
				steamIn.flowrate = sbr * bioFeedRate / steamIn.x[4];        // H2O is index 4 for Strean
			} else
			{
				// When steamIn.isMolar == true, both x and flowrate are molar
				// Calculate mass rate of steam from sbr, then convert to moles
				double wH2O = steamIn.x[4] * MW.H2O / steamIn.avgMW;
				double mDotTot = sbr * bioFeedRate / wH2O;
				steamIn.flowrate = mDotTot;		// Mass flow rate
				//steamIn.flowrate = mDotTot / steamIn.avgMW;
			}
		}

	}
}

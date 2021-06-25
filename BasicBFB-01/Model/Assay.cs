﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicBFB
{
	/* WtFrac Class
	 *		Used to manage values of mass fraction of individual elements
	 *		to model pyrolysis gas bulk composition.  Note the sum of mass fractions
	 *		for all elements must be equivalent to 1.0 by definition
	 *		
	 *		Fractions are stored in the w array with the following index to element mappings:
	 *			INDEX		ELEMENT
	 *			  0			 C	(Carbon)
	 *			  1			 H	(Hydrogen)
	 *			  2			 O	(Oxygen)
	 *			  3			 N	(Nitrogen)
	 *			  4			 S	(Sulfur)
	 *			  5			 X	(Inorganics/Ash components)
	 */

	public class Assay
	{
		public const int numElements = 6;
		public double[] w { get; set; }             // Elemental mass fractions central to this class

		public double p { get; set; }               // bar
		public double T { get; set; }               // °C
		public double flow { get; set; }		// Mass flow rate (might go somewhere else), kg/s


		// --------------------------------------------------------------------------
		// --- Next few params are used when working with PROXIMATE ANALYSIS results
		public double fracMoisture { get; set; }			// In weight fraction
		public double fracAsh { get; set; }				// In weight fraction
		public double fracVolatiles { get; set; }          // In weight fraction
		public double fracFixedCarbon { get; set; }        // In weight fraction
		public double feedParticleSize { get; set; }       // Measured in µm (avg diameter)
		// --------------------------------------------------------------------------
		

		// --------------------------------------------------------------------------
		// ----------------------------- CONSTRUCTORS -------------------------------
		public Assay()
		{
			this.w = new double[numElements];
			for (int i = 0; i < numElements; i++)
			{
				w[i] = 1.0 / numElements;
			}
		}

		// This constructor allows the programmer to initialize compositions
		public Assay(double[] wKnown)
		{
			// Initialize w so that it is of length numElements
			this.w = new double[numElements];
			int n = numElements;
			this.flow = 0.0;
			this.p = 1.01325;
			this.T = 25.0;

			fracMoisture = 0.0;
			fracAsh = 0.0;
			fracVolatiles = 0.0;
			fracFixedCarbon = 0.0;
			feedParticleSize = 0.0;

			// Check if wKnown length is less than  w.  If shorter, assume wKnown has
			// same index to element mappings
			// If wKnown is longer than numElements, w will contain the first numElements values from wKnown

			if (wKnown.Length <= numElements)
			{
				n = wKnown.Length;
			}

			for (int i = 0; i < n; i++)
			{
				w[i] = wKnown[i];
			}
		}


		// --------------------------------------------------------------------------
		// -------------------------------- METHODS ---------------------------------

		// This one makes a deep copy.  Since numElements is const,
		// we don't have to count elements in arrays
		public Assay Clone()
		{
			Assay wfNew = new Assay(this.w);
			wfNew.flow = this.flow;

			wfNew.fracMoisture = this.fracMoisture;
			wfNew.fracAsh = this.fracAsh;
			wfNew.fracVolatiles = this.fracFixedCarbon;
			wfNew.fracFixedCarbon = this.fracFixedCarbon;
			wfNew.feedParticleSize = this.feedParticleSize;
			return wfNew;
		}


		// Normalizes elemental fractions in this.w so they all sum to 1.0
		public void normalize()
		{
			double sum = 0.0;
			foreach (double v in w)
			{
				sum += v;
			}

			if (sum <= 0) {
				Console.WriteLine("WARNING: Could not normalize values in WtFrac - sum = {0}", sum);
				return; 
			}

			for (int i = 0; i < w.Length; i++)
			{
				w[i] = w[i] / sum;
			}
		}


		// Returns a new Assay object with elemental wt fractions normalized
		public Assay normalized()
		{
			Assay wfNormed = this.Clone();
			wfNormed.normalize();
			return wfNormed;
		}


		// Zeros out elements besides C, H, O and renomrmalize
		public void removeNonCHO()
		{
			for (int i = 3; i < numElements; i++)
			{
				w[i] = 0.0;
			}

			this.normalize();
		}


		// Without modifying this object, return the normalized array with CHO the only nonzero values
		// Length of return array is equal to numElements
		public Assay segregatedCHO()
		{
			Assay wfrac2 = new Assay(this.w);
			wfrac2.removeNonCHO();
			return wfrac2;
		}


		// Mathematically adjust wt fractions that result if all water driven off
		// Typically would normalize before running this step
		public void fullyDry()
		{
			// Since w[i] sums to 1, assume total mass of aliquot is 1.0 gram
			// This makes wtFrac values agree with pretend masses removed from "sample
			Assay dryWtFrac = this.Clone();
			double mH2O = this.fracMoisture;
			double mHRemoved = mH2O * 2.0 * MW.H / MW.H2O;
			double mORemoved = mH2O - mHRemoved;

			this.fracMoisture = 0.0;

			// Now subtract those evaporative mass losses from mass fraction array w
			this.w[1] -= mHRemoved;
			this.w[2] -= mORemoved;

			// Replace negatives with zero in case %moisture is overestimated
			if (this.w[1] < 0.0) { this.w[1] = 0.0; }	
			if (this.w[2] < 0.0) { this.w[2] = 0.0; }

			// Normalize both proximate and ultimate analysis values
			dryWtFrac.normalize();
			double wetProx = mH2O + this.fracAsh + this.fracVolatiles + this.fracFixedCarbon;
			double dryProx = wetProx - mH2O;
			double factor = wetProx / dryProx;
			this.fracAsh *= factor;
			this.fracVolatiles *= factor;
			this.fracFixedCarbon *= factor;
		}


		// This method is used to adjust biomass assay composition after 
		// remmoval of H2O, N2, O2 and S contaminants
		public Assay cleanedAndDried()
		{
			Assay wfClean = this.normalized();
			wfClean.fullyDry();
			wfClean.segregatedCHO();
			return wfClean;
		}
	}
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicBFB_01.Model
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

	class WtFrac
	{
		public const int numElements = 6;
		public double[] w { get; set; }

		public WtFrac()
		{
			this.w = new double[numElements];
			for (int i = 0; i < numElements; i++)
			{
				w[i] = 1.0 / numElements;
			}
		}


		public WtFrac(double[] wKnown)
		{
			// Initialize w so that it is of length numElements
			this.w = new double[numElements];
			int n = numElements;

			// Check if wKnown length is less than  w.  If shorter, assume wKnown has
			// same index to element mappings
			// If wKnown is longer than numElements, w will contain the first numElements values from wKnown

			if (wKnown.Length < numElements)
			{
				n = wKnown.Length;
			}

			for (int i = 0; i < n; i++)
			{
				w[i] = wKnown[i];
			}
			Console.WriteLine("Initialized a WtFrac object from input array - {0} elements copied.", n);
		}


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


		public double[] normalized()
		{
			double[] wNormed = new double[numElements];
			double sum = 0.0;
			foreach (double v in w)
			{
				sum += v;
			}

			if (sum <= 0)
			{
				Console.WriteLine("WARNING: Could not normalize values in WtFrac object - sum = {0}", sum);
				sum = 1.0;
			}

			for (int i = 0; i < numElements; i++)
			{
				wNormed[i] = this.w[i] / sum;
			}

			return wNormed;
		}


		// Modify the object's composition to consider only CHO in normalized wt%
		public void constrainToOnlyCHO()
		{
			// Zero out the non-CHO element wt fractions
			w[3] = 0.0;		// Nitrogen
			w[4] = 0.0;		// Sulfur
			w[5] = 0.0;     // Inorganics

			// Normalize the remaining values in w
			this.normalize();
		}


		// Without modifying this object, return the normalized array with CHO the only nonzero values
		// Length of return array is equal to numElements
		public WtFrac constrainedToOnlyCHO()
		{
			WtFrac wfrac2 = new WtFrac(this.w);
			wfrac2.constrainToOnlyCHO();
			return wfrac2;
		}
	}
}

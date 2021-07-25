using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Common = BasicBFB.Model.Common;

namespace BasicBFB.Model
{
	public class Effluent
	{
		public Stream total;			// All components in vapor (wet, mass basis)

		public Stream productGas;		// Only CO, CO2, CH4 and H2 (dry, mole basis)
		public double[] productYields;  // Wt frac. dry product gas, tar, and char
		public double mDotProducts = 0.0;

		public double mDotCharOut = 0.0;

		public double z = 0.0;			// Height in reactor where effluent is measured


		// --------------------------------------------------------------------------------
		//									CONSTRUCTORS
		// --------------------------------------------------------------------------------

		// Default constructor
		public Effluent()
		{
			double xi = 1.0 / ((double)Stream.numComp);
			double[] x0 = new double[Stream.numComp];
			for (int i = 0; i < Stream.numComp; i++)
			{
				x0[i] = xi;
			}

			this.total = new Stream(x0);

			xi = 0.25;
			for (int i = 0; i < 4; i++)
			{
				x0[i] = xi;
			}

			for (int i = 4; i < Stream.numComp; i++)
			{
				x0[i] = 0.0;
			}

			this.productGas = new Stream(x0, isMolar: true);
			this.productYields = new double[3];
		}


		// Constructor from mass rates in total vapor
		public Effluent(double[] mDot, double zOut, GasifierParams param, double mDotCharOut = 0.0)
		{
			this.z = zOut;
			this.mDotCharOut = mDotCharOut;
			this.total = Stream.CreateFromMassRates(mDot, param.T, param.p, false);

			double[] mDotProductGas = new double[Stream.numComp];
			Array.Copy(mDot, mDotProductGas, 4);
			this.productGas = Stream.CreateFromMassRates(mDotProductGas, param.T, param.p, true);

			this.productYields = new double[3];
			this.mDotProducts = productGas.mDotTotal + mDot[(int)Common.Component.Tar] + mDotCharOut;
			productYields[(int)Common.Phase.Gas] = productGas.mDotTotal / mDotProducts;
			productYields[(int)Common.Phase.Tar] = mDot[(int)Common.Component.Tar] / mDotProducts;
			productYields[(int)Common.Phase.Char] = mDotCharOut / mDotProducts;
		}

	}
}

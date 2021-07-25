using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using BasicBFB.Model.Common;


namespace BasicBFB.Model
{
	public class ReactorFreeboard
	{
		public GasifierParams param;
		public Effluent bedEffluentIn { get; init; }
		
		// Final reactor effluent
		public Effluent effluent { get; private set; }

		private double[] zSpan { get; init; }     // {Lbed, Ltot}
		private double[] mz0 { get; init; }       // Initial condition = Mass rates from reactor bed out

		public List<double> zList { get; private set; }
		public List<double[]> mDotList { get; private set; }
		// j = CO, CO2, CH4, H2, H2O, N2, O2, S, Tar


		// --------------------------------------------------------------------------------
		//									CONSTRUCTOR
		// --------------------------------------------------------------------------------
		public ReactorFreeboard(GasifierParams param, Effluent bedEffluent)
		{
			this.param = param;
			this.bedEffluentIn = bedEffluent;
			this.effluent = new Effluent();
			this.zSpan = new double[2] { bedEffluentIn.z, param.Ltot };
			this.mz0 = bedEffluent.total.mDot;
			this.zList = new List<double>(200);
			this.mDotList = new List<double[]>(200);
		}


		// --------------------------------------------------------------------------------
		//							  GENERAL HELPER FUNCTIONS
		// --------------------------------------------------------------------------------

		// Gas phase mole fraction for specified component
		// at z position corresponding to specified zIndex
		private double yMolar(Component comp, in double[] mz)
		{
			double sumMoles = 0.0;

			for (int j = 0; j < Stream.numComp; j++)
			{
				sumMoles += mz[j] / MW.all[j];
			}

			double y = mz[(int)comp] / MW.all[(int)comp];
			y /= sumMoles;

			return y;
		}

		// Partial pressure of specified component in gas phase 
		// at z position for specified index (bara)
		private double partialP(Component comp, in double[] mz)
		{
			return param.p * yMolar(comp, mz);
		}

		// Concentration of specified component in gas at zIndex [=] mol/m3
		private double conc(Component comp, in double[] mz)
		{
			double Tkelvin = param.T + 273.15;
			double c = 1.0e5 * partialP(comp, mz) / (Const.Rgas * Tkelvin);
			return c;
		}

		// Superficial gas phase velocity at z[zIndex]
		private double gasU(double[] mz)
		{
			double Tkelvin = param.T + 273.15;
			double u = 0.0;

			for (int j = 0; j < Stream.numComp; j++)
			{
				u += mz[j] / MW.all[j];
			}

			u *= (Const.Rgas * Tkelvin) / (1.0e5 * param.p * param.Axs);
			return u;
		}


	}
}

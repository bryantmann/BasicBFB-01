using BasicBFB.Model.Common;

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicBFB.Model
{
	public class ReactorZone
	{
		public GasifierParams param;
		public Pyrolysis pyro;

		protected double[] zSpan { get; init; }		// Range of z values for use by ODE solver	
		protected double[] mz0 { get; set; }        // Initial condition = Mass rates from reactor bed out

		public List<double> zList;						// Heights z[i] mDotList[i] 
		public List<double[]> mDotList;					// Mass flow rate solutions to ODE
		public Effluent effluent { get; protected set; }            // Section's effluent stream data


		// --------------------------------------------------------------------------------
		//					MAIN CALCULATION ALGORITHMS (to be overriden)
		// --------------------------------------------------------------------------------

		public virtual Effluent calcEffluent()
		{
			solve();
			return effluent;
		}
		public virtual void solve()
		{
			Console.WriteLine("Called solve() method from ReactorZone framework class!");
			Console.WriteLine("Ensure the child classes correctly override the solve() method");
		}

		protected virtual double[] dmDot(double z, in double[] mz)
		{
			double[] nothing = new double[0];
			Console.WriteLine("Called dmDot() method from ReactorZone framework class!");
			Console.WriteLine("Ensure the child classes correctly override the dmDot() method");
			return nothing;
		}


		// --------------------------------------------------------------------------------
		//							  GENERAL HELPER FUNCTIONS
		// --------------------------------------------------------------------------------

		// Gas phase mole fraction for specified component
		// at z position corresponding to specified zIndex
		protected double yMolar(Component comp, in double[] mz)
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
		protected double partialP(Component comp, in double[] mz)
		{
			return param.p * yMolar(comp, mz);
		}

		// Concentration of specified component in gas at zIndex [=] mol/m3
		protected double conc(Component comp, in double[] mz)
		{
			double Tkelvin = param.T + 273.15;
			double c = 1.0e5 * partialP(comp, mz) / (Const.Rgas * Tkelvin);
			return c;
		}

		// Superficial gas phase velocity at z[zIndex]
		public double gasU(in double[] mz)
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

		// Local bed voidage eps(z) based on zIndex
		internal virtual double epsilon(double z, in double[] mz)
		{
			Console.WriteLine("Called epsilon() method from ReactorZone framework class!");
			Console.WriteLine("Ensure the child classes correctly override the epsilon() method");
			return 0.0;
		}


		// --------------------------------------------------------------------------------
		//						REACTION RATE HELPERS - To be overriden
		// --------------------------------------------------------------------------------

		// Gasification reaction rates returned in array
		// Order is rC1, rC2, rC3, rSMR, rWGS
		// Indexed the same as Reaction enumeration raw values
		internal virtual double[] gasifyRates(double z, in double[] mz)
		{
			double[] nothing = new double[0];
			Console.WriteLine("Called gasifyRates() method from ReactorZone framework class!");
			Console.WriteLine("Ensure the child classes correctly override the gasifyRates() method");
			return nothing;
		}

		// CO rate expression (dmCO/dz)
		internal virtual double dmCO(in double[] mz, in double[] gasRates, double U, double eps)
		{
			Console.WriteLine("Called dmCO() method from ReactorZone framework class!");
			Console.WriteLine("Ensure the child classes correctly override the dmCO() method");
			return 0.0;
		}

		// CO2 rate expression (dmCO2/dz)
		internal virtual double dmCO2(in double[] mz, in double[] gasRates, double U, double eps)
		{
			Console.WriteLine("Called dmCO2() method from ReactorZone framework class!");
			Console.WriteLine("Ensure the child classes correctly override the dmCO2() method");
			return 0.0;
		}

		// CH4 rate expression (dmCH4/dz)
		internal virtual double dmCH4(in double[] mz, in double[] gasRates, double U, double eps)
		{
			Console.WriteLine("Called dmCH4() method from ReactorZone framework class!");
			Console.WriteLine("Ensure the child classes correctly override the dmCH4() method");
			return 0.0;
		}

		// H2 rate expression (dmH2/dz)
		internal virtual double dmH2(in double[] mz, in double[] gasRates, double U, double eps)
		{
			Console.WriteLine("Called dmH2() method from ReactorZone framework class!");
			Console.WriteLine("Ensure the child classes correctly override the dmH2() method");
			return 0.0;
		}

		// H2O rate expression (dmH2O/dz)
		internal virtual double dmH2O(in double[] mz, in double[] gasRates, double U, double eps)
		{
			Console.WriteLine("Called dmH2O() method from ReactorZone framework class!");
			Console.WriteLine("Ensure the child classes correctly override the dmH2O() method");
			return 0.0;
		}

		// Tar rate expression (dmTar/dz)
		internal virtual double dmTar(in double[] mz, double U, double eps)
		{
			Console.WriteLine("Called dmTar() method from ReactorZone framework class!");
			Console.WriteLine("Ensure the child classes correctly override the dmTar() method");
			return 0.0;
		}

		// Char rate expression (dmChar/dz)
		internal virtual double dmChar(in double[] mz, double[] gasRates, double U)
		{
			Console.WriteLine("Called dmChar() method from ReactorZone framework class!");
			Console.WriteLine("Ensure the child classes correctly override the dmChar() method");
			return 0.0;
		}

	}
}

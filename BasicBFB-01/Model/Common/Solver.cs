using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicBFB.Model.Common
{
	public delegate double Del1(double x);
	
	
	public static class Solver
	{
		// --------------------------------------------------------------------------------
		// Parameters for the Dormand-Prince numerical ODE algorithm (used in MATLAB's ode45)
		// --------------------------------------------------------------------------------

		//public static double[,] B2 = { {0.2, (3.0 / 40.0), (44.0/45.0), (19372.0 / 6561.0), (9017.0 / 3168.0), (35.0 / 384.0)},
		//							  {0.0, (9.0 / 40.0), (-56.0 / 15.0), (-25360.0 / 2187.0), (-355.0 / 33.0), 0.0},
		//							  {0.0, 0.0, (32.0 / 9.0), (64448.0 / 6561.0), (46732.0 / 5247.0), (500.0 / 1113.0)},
		//							  {0.0, 0.0, 0.0, (-212.0 / 729.0), (49.0 / 176.0), (125.0 / 192.0)},
		//							  {0.0, 0.0, 0.0, 0.0, (-5103.0 / 18656.0), (-2187.0 / 6784.0)},
		//							  {0.0, 0.0, 0.0, 0.0, 0.0, (11.0 / 84.0)},
		//							  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
		//							};

		// --------------------------------------------------------------------------------



		// Ridder's root finding method
		// Translated from Numerical Recipes in C
		public static double ridder(Del1 func, double[] xbounds, double xacc)
		{
			const int MAXIT = 60;
			const double UNUSED = (-1.11e30);

			double ans, fh, fl, fm, fnew, s, xh, xl, xm, xnew;

			fl = func(xbounds[0]);
			fh = func(xbounds[1]);

			if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0))
			{
				xl = xbounds[0];
				xh = xbounds[1];
				ans = UNUSED;

				for (int j = 0; j < MAXIT; j++)
				{
					xm = 0.5 * (xl - xh);
					fm = func(xm);
					s = Math.Sqrt(fm * fm - fl * fh);
					if (s == 0.0) { return ans; }

					xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);

					if (Math.Abs(xnew - ans) <= xacc) { return ans; }

					ans = xnew;
					fnew = func(ans);

					if (fnew == 0.0) { return ans; }

					if (sign(fm, fnew) != fm)
					{
						xl = xm;
						fl = fm;
						xh = ans;
						fh = fnew;
					}
					else if (sign(fl, fnew) != fl)
					{
						xh = ans;
						fh = fnew;
					}
					else if (sign(fh, fnew) != fh)
					{
						xl = ans;
						fl = fnew;
					}
					else
					{
						throw new NotConvergedException(xl, xh);
					}
				}
			}
			else
			{
				if (fl == 0.0) { return xbounds[0]; }
				if (fh == 0.0) { return xbounds[1]; }

				// Getting here means xbounds[0] and [1] are nonzero with the same sign
				string message = "Root must be bracketed in Solver.ridder()";
				string argname = "xbounds";
				throw new ArgumentException(message, argname);
			}

			return 0.0;			// Cannot get here 
		}


		public static double sign(double a, double b)
		{
			double ans = Math.Abs(a);
			if (b < 0.0) { ans = -ans; }
			return ans;
		}

	}
}

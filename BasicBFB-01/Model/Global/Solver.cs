using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicBFB
{
	public delegate double Del1(double x);	
	
	
	public static class Solver
	{
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

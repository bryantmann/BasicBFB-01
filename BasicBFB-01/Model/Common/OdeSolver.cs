using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicBFB.Model.Common
{
	public delegate double[] OdeDel(double t, double[] y);

	public class OdeSolver
	{
		// --------------------------------------------------------------------------------
		// Parameters for the Dormand-Prince numerical ODE algorithm (used in MATLAB's ode45)
		// --------------------------------------------------------------------------------
		public const double pow = 0.2;

		public static double[] A = { 0.0, 0.2, 0.3, 0.8, (8.0 / 9.0), 1.0, 1.0 };

		public static double[,] B = {   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
										{ (1.0 / 5.0), 0.0, 0.0, 0.0, 0.0, 0.0 },
										{ (3.0 / 40.0), (9.0 / 40.0), 0.0, 0.0, 0.0, 0.0 },
										{ (44.0 / 45.0), (-56.0 / 15.0), (32.0 / 9.0), 0.0, 0.0, 0.0 },
										{ (19372.0 / 6561.0), (-25360.0 / 2187.0), (64448.0 / 6561.0), (-212.0 / 729.0), 0.0, 0.0 },
										{ (9017.0 / 3168.0), (-355.0 / 33.0), (46732.0 / 5247.0), (49.0 / 176.0), (-5103.0 / 18656.0), 0.0 },
										{ (35.0 / 384.0), 0.0, (500.0 / 1113.0), (125.0 / 192.0), (-2187.0 / 6784.0), (11.0 / 84.0) }
									};

		public static double[,] C = {   { (35.0 / 384.0), 0.0, (500.0 / 1113.0), (125.0 / 192.0), (-2187.0 / 6784.0), (11.0 / 84.0), 0.0 },
										{ (5179.0 / 57600.0), 0.0, (7571.0 / 16695.0), (393.0 / 640.0), (-92097.0 / 339200.0), (187.0 / 2100.0), (1.0 / 40.0) }
									};

		public static double[] E = { (71.0 / 57600.0), 0.0, (-71.0 / 16695.0), (71.0 / 1920.0),
									 (-17253.0 / 339200.0), (22.0 / 525.0), (-1.0 / 40.0) };
		// --------------------------------------------------------------------------------
		
		// Solver parameters
		public double absTol = 1.0e-6;
		public double relTol = 1.0e-3;
		public double hMaxFrac = 0.1;
		public double absh = 0.01;
		public double tDir = 1.0;
		public double h
		{
			get
			{
				return tDir * absh;
			}
		}


		// Algorithm based on MATLAB's ode45
		public (List<double> t, List<double[]> Y) rk45(OdeDel dy, double[] tspan, double[] y0)
		{
			double deltaT = Math.Abs(tspan[1] - tspan[0]);
			tDir = (tspan[1] > tspan[0]) ? 1.0 : -1.0;

			List<double> tOut = new List<double>(200);
			List<double[]> Y5 = new List<double[]>(200);

			double t = tspan[0];
			double[] y = y0;

			double hmin = 16.0 * getEps(t);
			double hmax = hMaxFrac * deltaT;

			// Calculate an initial step size using dy(y, t)
			absh = Math.Min(0.01 * deltaT, hmax);

			double[] f0 = dy(t, y0);
			double rh = getInitialStepSize(y, f0);
			if (absh*rh > 1.0)
			{
				absh = 1.0 / rh;
			}
			absh = Math.Max(absh, hmin);

			// Initialize first row of Y and first element of t
			tOut.Add(t);
			Y5.Add(y);

			double tNext = 0.0;
			double err = 0.0;
			double[][] V = new double[6][];
			double[] y4Next = new double[y.Length];
			double[] yi = new double[y.Length];
			double[] yErr = new double[y.Length];
			double[] maxAbsY = new double[y.Length];
			double[] yErrRatio = new double[y.Length];
			double ti = 0.0;
			double factor = 0.0;
			double temp = 0.0;

			bool isDone = false;
			while (!isDone)
			{
				int n = tOut.Count - 1;      // index of last element
				t = tOut[n];
				y = Y5[n];

				//double[] y4Next = new double[y.Length];
				double[] y5Next = new double[y.Length];

				// Get step size, constrained by hmin based on last t value
				hmin = 16.0 * getEps(t);
				absh = Math.Min(hmax, Math.Max(hmin, absh));

				// Stretch the step if within 10% of tfinal - t
				if (1.1 * absh >= Math.Abs(tspan[1] - t))
				{
					absh = Math.Abs(tspan[1] - t);
					isDone = true;
				}

				// -- Loop for advancing one step, continuing only if error is within tolerance
				bool notFailed = true;
				while (true)
				{
					// Calculate the six intermediate slopes 
					
					V[0] = dy(t, y);
					
					for (int i = 1; i <= 5; i++)
					{
						ti = t + A[i] * h;
						
						for (int k = 0; k < y.Length; k++)
						{
							double sum_k = 0.0;
							for (int j = 0; j < i; j++)
							{
								sum_k += h * B[i, j] * V[j][k];
							}
							yi[k] = y[k] + sum_k;
						}

						V[i] = dy(ti, yi);
					}

					// Update value for tNext using current value of h
					tNext = t + h;

					// Using the slopes, calculate y(n+1) values using both 
					// the 4th and 5th order formulae
					for (int k = 0; k < y.Length; k++)
					{
						// 4th order formula
						double sum4_k = 0.0;
						for (int i = 0; i <= 4; i++)
						{
							sum4_k += h * C[0, i] * V[i][k];
						}

						// 5th order formula
						double sum5_k = 0.0;
						for (int i = 0; i <= 5; i++)
						{
							sum5_k += h * C[1, i] * V[i][k];
						}

						y4Next[k] = y[k] + sum4_k;
						y5Next[k] = y[k] + sum5_k;
					}

					// Estimate the error with current h using y4Next and y5Next values
					for (int k = 0; k < y.Length; k++)
					{
						yErr[k] = y5Next[k] - y4Next[k];
					}

					err = calcError(y, y5Next, yErr, ref maxAbsY, ref yErrRatio);

					// Accept this solution only if weighted error err is no more	
					// than the tolerance relTol.  Estimate an h that will yield
					// an error of relTol on the next step/try, and use 0.8 of this
					// value to avoid failures
					if (err > relTol)
					{
						if (absh <= hmin)
						{
							string message = $"OdeSolver.rk45 error at t = {t:g4} exceeds tolerance even with minimum step size of {hmin:g4}.";
							throw new NotConvergedException(t, err, message);
						}

						if (notFailed)
						{
							notFailed = false;
							//double factor = 0.8 * Math.Pow((relTol / err), pow);
							factor = 0.8 * Math.Pow((relTol / err), pow);
							factor = Math.Max(0.1, factor);
							absh = Math.Max(hmin, absh * factor);
						} else
						{
							absh = Math.Max(hmin, 0.5 * absh);
						}

						isDone = false;
					} else
					{
						break;
					}
				}	// Close of inner while loop for single step



				if (isDone)
				{
					break;
				}

				// If there were no failures compute a new h
				if (notFailed)
				{
					// Note that absh may shrink by 0.8, and that err may be 0.
					//double temp = 1.25 * Math.Pow((err / relTol), pow);
					temp = 1.25 * Math.Pow((err / relTol), pow);
					if (temp > 0.2)
					{
						absh = absh / temp;
					} else
					{
						absh = 5.0 * absh;
					}
				}

				// Append new values for tNext and y5Next to tOut and Y5
				tOut.Add(tNext);
				Y5.Add(y5Next);
			}

			//tOut.TrimExcess();
			//Y5.TrimExcess();

			return (tOut, Y5);
		}


		private double getInitialStepSize(in double[] y0, in double[] f0)
		{
			double threshold = absTol / relTol;
			double[] v = new double[y0.Length];
			for (int i = 0; i < y0.Length; i++)
			{
				v[i] = f0[i] / Math.Max(Math.Abs(y0[i]), threshold);
			}

			return normi(v) / (0.8 * Math.Pow(relTol, pow));
		}

		private static double norm(in double[] v)
		{
			double sumSq = 0.0;
			for (int i = 0; i < v.Length; i++)
			{
				sumSq += v[i] * v[i];
			}

			return Math.Sqrt(sumSq);
		}

		private static double normi(in double[] v)
		{
			double maxabs = Math.Abs(v[0]);
			for (int i = 1; i <  v.Length; i++)
			{
				maxabs = Math.Max(maxabs, Math.Abs(v[i]));
			}

			return maxabs;
		}

		private double getEps(double x)
		{
			if (x == 0.0)
			{
				x += Double.Epsilon;
			}
			double machEps = x;
			do
			{
				machEps /= 10.0d;
			}
			while ((double)(x + machEps) != x);
			return machEps;
		}


		private double calcError(in double[] yOld, in double[] yNew, in double[] yE, ref double[] maxAbsY, ref double[] yErrRatio)
		{
			double threshold = absTol / relTol;
			double err = 0.0;

			//maxAbsY = new double[yOld.Length];
			for (int i = 0; i < yNew.Length; i++)
			{
				double maxYi = Math.Max(Math.Abs(yOld[i]), Math.Abs(yNew[i]));
				maxAbsY[i] = Math.Max(maxYi, threshold);
			}

			//yErrRatio = new double[yOld.Length];
			for (int i = 0; i < yOld.Length; i++)
			{
				yErrRatio[i] = yE[i] / maxAbsY[i];
			}

			err = normi(yErrRatio);
			return err;
		}
	}
}

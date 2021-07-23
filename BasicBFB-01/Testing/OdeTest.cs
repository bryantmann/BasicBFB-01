using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

using BasicBFB.Model.Common;

namespace BasicBFB.Testing
{
	public class OdeTest
	{
		public double[] tRange = { 0.0, 30.0 };

		public double[] springODE(double t, double[] y)
		{
			double[] dy = new double[2];
			dy[0] = y[1];		// dx/dt = v
			dy[1] = -y[0];		// dv/dt = -x
			return dy;
		}

		public double[] springAnal(double t)
		{
			double[] y = new double[2];
			y[0] = Math.Cos(t);
			y[1] = -Math.Sin(t);
			return y;
		}


		public void springTest()
		{
			string dash20 = new String('=', 20);
			Console.WriteLine(dash20 + ">  ODE TEST: Spring  <" + dash20);

			// 1. Numerical solution using rk45
			Console.Write("\n--->  Numerical Solution using rk45...  ");
			double[] y0 = { 1.0, 0.0 };
			OdeDel ode = springODE;
			OdeSolver odeSolver = new OdeSolver();
			(List<double> t1, List<double[]> Y1) = odeSolver.rk45(ode, tRange, y0);
			string csv1 = makeCSV(t1, Y1);
			Console.WriteLine("Complete!");

			// 2. Analytical solution
			Console.Write("\n--->  Analytical Solution...  ");
			const double dt2 = 0.05;
			int N2 = (int)(Math.Round( Math.Abs(tRange[1] - tRange[0]) / dt2 ) );
			List<double> t2 = new List<double>();
			List<double[]> Y2 = new List<double[]>();
			for (int i = 0; i <= N2; i++)
			{
				double ti = i * dt2;
				double[] yi = springAnal(ti);
				t2.Add(ti);
				Y2.Add(yi);
			}
			string csv2 = makeCSV(t2, Y2);
			Console.WriteLine("Complete!");

			// Save results in two files in my temp directory
			Console.Write("\n--->  Writing results to files in ~/temp/...  ");
			writeStrToFile(csv1, "springNumerical.csv");
			writeStrToFile(csv2, "springAnalytical.csv");
			Console.WriteLine("Complete!");
		}


		public string makeCSV(List<double> t, List<double[]> Y)
		{
			int n = t.Count;
			string s = "";

			// header stuff goes here
			string hdr = "t,x,v\n";
			s += hdr;

			for (int i = 0; i < n; i++)
			{
				double tObs = t[i];
				double[] yObs = Y[i];
				string line = $"{tObs},{yObs[0]},{yObs[1]}\n";
				s += line;
			}
			return s;
		}


		public static void writeStrToFile(string str, string filename)
		{
			string path = $"C:\\Users\\bryan\\temp\\{filename}";
			File.WriteAllText(path, str);
		}


		public void printCSV(List<double> t, List<double[]> Y)
		{
			int n = t.Count;
			string s = "";

			// header stuff goes here
			string dashes = new String('-', 80);
			string hdr = "t,x,v\n";
			s += dashes + "\n";
			s += hdr;

			for (int i = 0; i < n; i++)
			{
				double tObs = t[i];
				double[] yObs = Y[i];
				string line = $"{tObs:F3},{yObs[0]:F3},{yObs[1]:F3}\n";
				s += line;
			}

			s += dashes + "\n";

			Console.WriteLine(s);
		}

	}
}

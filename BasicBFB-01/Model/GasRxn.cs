using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using BasicBFB.Model.Common;

namespace BasicBFB.Model
{
	static public class GasRxn
	{
		// Reaction enthalpies, all with units of kJ / mol
		// Use with Reaction enum to get index for a reaction type
		static double[] dHrxn = { 172.0, 131.0, -75.0, 206.0, -41.0 };

		// Input parameter units:
		//		Tkelvin [=]  K
		//		pi [=] bara
		//		Ci [=] mol/m3

		// Rate of Boudouard reaction [s^-1]
		//		C(s) + CO2  -->  2CO		dH° = +172 kJ/mol
		static public double rateC1(double pCO2, double Tkelvin)
		{
			//double Tkelvin = T + 273.15;
			double r = 3.1e6 * Math.Exp(-215000.0 / (Const.Rgas * Tkelvin));
			r *= Math.Pow(pCO2, 0.38);
			return r;
		}

		// Rate of water-gas reaction [s^-1]
		//		C(s) + H2O  -->  CO + H2	dH° = +131 kJ/mol
		static public double rateC2(double pH2O, double Tkelvin)
		{
			//double Tkelvin = T + 273.15;
			double r = 2.62e8 * Math.Exp(-237000.0 / (Const.Rgas * Tkelvin));
			r *= Math.Pow(pH2O, 0.57);
			return r;
		}

		// Rate of methanation reaction [s^-1]
		//		C(s) + 2H2  -->  CH4		dH° = -75 kJ/mol
		static public double rateC3(double pH2, double Tkelvin)
		{
			//double Tkelvin = T + 273.15;
			double pH2_MPa = pH2 / 10.0;
			double r = 16.4 * Math.Exp(-94800.0 / (Const.Rgas * Tkelvin));
			r *= Math.Pow(pH2_MPa, 0.93);
			return r;
		}

		// Rate of steam methane reforming (SMR) reaction [mol/m3.s]
		//		CH4 + H2O  <==>  CO + 3H2 
		static public double rateSMR(double cCH4, double cH2O, double Tkelvin)
		{
			//double Tkelvin = T + 273.15;
			double r = 3.0e5 * Math.Exp(-125000.0 / (Const.Rgas * Tkelvin));
			r *= cCH4 * cH2O;
			return r;
		}

		// Rate of water gas shift reaction (WGS) [mol/m3.s]
		//		CO + H2O  <==>  CO2 + H2
		static public double rateWGS(double cCO, double cCO2, double cH2, double cH2O, 
										double Tkelvin)
		{
			//double Tkelvin = T + 273.15;
			double Keq = 0.0265 * Math.Exp(3968.0 / Tkelvin);

			double cTerm = cCO * cH2O - cCO2 * cH2 / Keq;

			double r = 2.78 * Math.Exp(-1510.0 / Tkelvin) * cTerm;
			return r;
		}

		//static public Dictionary<string, double> dHrxn2 = new Dictionary<string, double>
		//{
		//	{ "C1",  172.0 }, { "C2",  131.0 }, { "C3", -75.0 },
		//	{ "SMR", 206.0 }, { "WGS", -41.0 }
		//};
	}
}

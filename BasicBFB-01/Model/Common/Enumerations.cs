using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicBFB.Model.Common
{
	public enum Phase: int
	{
		Gas = 0,
		Tar = 1,
		Char = 2
	}


	// Raw values coincide with element index in Assay's w array (mass fractions)
	// Index order is used consistently throughout model (well, it's supposed to be)
	public enum Element: int
	{
		C = 0,
		H = 1,
		O = 2,
		N = 3,
		S = 4
	}


	// Raw values coincide with component index in Stream objects and elsewhere
	public enum Component: int
	{
		CO = 0,
		CO2 = 1,
		CH4 = 2,
		H2 = 3,
		H2O = 4,
		N2 = 5,
		O2 = 6,
		S = 7,
		Tar = 8
	}


	// Gasification reaction types - Raw values coincide with array index
	public enum Reaction: int
	{
		C1 = 0,
		C2 = 1, 
		C3 = 2, 
		SMR = 3,
		WGS = 4
	}
}

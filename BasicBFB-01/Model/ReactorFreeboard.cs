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

	}
}

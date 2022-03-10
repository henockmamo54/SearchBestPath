// Labeling_Path.h

#pragma once

using namespace System;
using namespace System::Collections::Generic;

namespace Labeling_Path {

	public ref class Label_Path
	{
		// TODO: Add your methods for this class here.

	public:
		array <float>^ fFractonalSynthesis, ^ fBWE_Array;
		float fBWE;

		List<String^>^ Labeling_Path_Fractional_Synthesis(array <float, 2>^ Fractional_Synthesis_Rate,
			int Nraw, int Ncol, int NEH, float fBWE);


		List<String^>^ Labeling_Path_Fractional_Synthesis2(array <float, 2>^ Fractional_Synthesis_Rate,
			int Nraw, int Ncol, int NEH, float fBWE);


		float Diff_Fraction_Rates(float Fractional_Synthesis_Rate1,
			float Fractional_Synthesis_Rate2);


	};
}

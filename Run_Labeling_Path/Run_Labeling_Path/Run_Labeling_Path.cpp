// Run_Labeling_Path.cpp : main project file.

#include "stdafx.h"
#include <stdlib.h>  
#include <stdio.h> 

using namespace System;
using namespace Labeling_Path;
using namespace System::Collections::Generic;

int main(array<System::String^>^ args)
{
	int i, j;

	float range_max, range_min, random_number;

	int NLim = 7, NProt = 20;;

	array <float, 2>^ mock_Time_Course = gcnew array <float, 2>(NProt, NLim);

	array <float, 1>^ time = gcnew array <float, 1>(NLim);

	range_max = .5; range_min = -.5;

	Label_Path^ getPath = gcnew Label_Path();

	for (i = 0; i < NProt; i++)
	{
		for (j = 0; j < NLim; j++)
		{
			random_number = (float)rand() / (RAND_MAX + 1) * (range_max - range_min)
				+ range_min;

			mock_Time_Course[i, j] = 1. - Math::Exp(-0.1 * j * 2) + random_number;
		}
	}

	for (i = 0; i < NProt; i++)
	{
		for (j = 0; j < NLim; j++)
		{
			printf("mock_Time[%d, %d] = %10.5f ", i, j, mock_Time_Course[i, j]);

		}
		printf("\n");
	}

	List<String^>^ values = getPath->Labeling_Path_Fractional_Synthesis(mock_Time_Course, NProt, NLim, 0, 0.);

	return 0;
}

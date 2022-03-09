// This is the main DLL file.

#include "stdafx.h"
#include "stdio.h"
#include "math.h"
#include "string"
#include "Labeling_Path.h"

using namespace System;
using namespace Labeling_Path;
using namespace System::Collections::Generic;

/*
*
*   The method computes the optimal labeling path given the
*   the time course of the fractional synthesis rate from
*   several proteins of a peptide. The path means that the
*   program chooses the best fractional synthesis rate time
*   points so that the overall sequence of the fractional
*   synthesis rate is monotonically increasing.
*   The input data is in Fractional_Synthesis_Rate
*   It is assumed that the columns are peptides,
*   and the rows are the time course values.
*
*/
List<String^>^ Label_Path::Labeling_Path_Fractional_Synthesis(array <float, 2>^ Fractional_Synthesis_Rate,
	int Nrow, int Ncol, int NEH, float fBWE)
{
	int i, j, i1, j1, i0, j0, iMax;

	array <float, 2>^ Score_Fractional_Synthesis_Rate_Time = gcnew array <float, 2>(Nrow, Ncol);

	array <float, 2>^ fBest_Fractional_Synthesis_Rate_Time = gcnew array <float, 2>(Nrow, Ncol);

	array <int, 2>^ index_X = gcnew array <int, 2>(Nrow, Ncol);

	array <int, 2>^ index_Y = gcnew array <int, 2>(Nrow, Ncol);

	List < String^ >^ bestpaths = gcnew List<String^>();

	List<float>^ bestpaths_score = gcnew List<float>();
	List<String^>^ bestpaths_and_score = gcnew List<String^>();

	float fSum, ftemp, fMinSlope;

	for (i = 0; i < Nrow; i++)
		for (j = 0; j < Ncol; j++)
			Score_Fractional_Synthesis_Rate_Time[i, j] = 0; // -1. * (Ncol + Nrow);

	for (i = 0; i < Ncol; i++)
	{
		//index_X[i, 0] = i;
		index_X[0, i] = i;
	}

	for (i = 0; i < Nrow; i++)
	{
		fBest_Fractional_Synthesis_Rate_Time[i, 0] = Fractional_Synthesis_Rate[i, 0];
	}

	for (i = 0; i < Nrow; i++)
		//for (i = 0; i < 1; i++)
	{
		for (j = 1; j < Ncol; j++)
		{
			fMinSlope = 100.0;
			//
			//frist checks to see if the next point from this peptide's time course is monotonic
			//if yes, then it will accept this value for the peptides' time course.
			//
			//if (fBest_Fractional_Synthesis_Rate_Time[i, j - 1] == Fractional_Synthesis_Rate[i, j - 1] &&
				//Fractional_Synthesis_Rate[i, j] > Fractional_Synthesis_Rate[i, j - 1] )

			////if (Fractional_Synthesis_Rate[i, j] <= 0) {
			////	fBest_Fractional_Synthesis_Rate_Time[i, j] = NAN;
			////	index_X[i, j] = -1;
			////	index_Y[i, j] = -1;
			////	continue;
			////}

			////// to hold the previous none zero fBest_Fractional_Synthesis_Rate_Time
			////float fBest_Fractional_Synthesis_Rate_Time_previous = fBest_Fractional_Synthesis_Rate_Time[i, j - 1];
			////float counter = 1;
			////while (fBest_Fractional_Synthesis_Rate_Time_previous == NAN)
			////{
			////	fBest_Fractional_Synthesis_Rate_Time_previous = fBest_Fractional_Synthesis_Rate_Time[i, j - 1 - counter];
			////	counter++;
			////}


			if (Fractional_Synthesis_Rate[i, j] > fBest_Fractional_Synthesis_Rate_Time[i, j - 1])

			{

				Score_Fractional_Synthesis_Rate_Time[i, j] = Score_Fractional_Synthesis_Rate_Time[i, j - 1] +
					Diff_Fraction_Rates(Fractional_Synthesis_Rate[i, j - 1], Fractional_Synthesis_Rate[i, j]);

				index_X[i, j] = i;

				index_Y[i, j] = j;

				fBest_Fractional_Synthesis_Rate_Time[i, j] = Fractional_Synthesis_Rate[i, j];
			}
			else
			{
				for (i0 = 0; i0 < Nrow; i0++)
				{
					//fSum = Score_Fractional_Synthesis_Rate_Time[i0, j - 1] +
						//Diff_Fraction_Rates(Fractional_Synthesis_Rate[i, j - 1], Fractional_Synthesis_Rate[i0, j],
							//NEH, fBWE);

					//fSum = Score_Fractional_Synthesis_Rate_Time[i0, j - 1] + 
					fSum = Score_Fractional_Synthesis_Rate_Time[i, j - 1] +
						//Diff_Fraction_Rates(Fractional_Synthesis_Rate[index_X[i, j - 1], index_Y[i, j - 1]], Fractional_Synthesis_Rate[i0, j]);

						Diff_Fraction_Rates(fBest_Fractional_Synthesis_Rate_Time[i, j - 1], Fractional_Synthesis_Rate[i0, j]);

					//if(7 == i)
					//printf("%d %d, %3.2f %3.2f %3.2f %3.2f %4.3f %3.2f;  ", i0, j, fBest_Fractional_Synthesis_Rate_Time[i, j - 1], Fractional_Synthesis_Rate[i0, j], 
						//fSum, Score_Fractional_Synthesis_Rate_Time[i, j-1],
						//Diff_Fraction_Rates(fBest_Fractional_Synthesis_Rate_Time[i, j - 1], Fractional_Synthesis_Rate[i0, j] ),
						//Score_Fractional_Synthesis_Rate_Time[i, j]);
					if (fSum >= Score_Fractional_Synthesis_Rate_Time[i, j])
					{

						//when multiple paths are possible, choose the path with the least slope
						// this should always find a solution, if it exists.
						//if (fMinSlope > Math::Abs(Fractional_Synthesis_Rate[i0, j] - fBest_Fractional_Synthesis_Rate_Time[i, j - 1]))
						//if(Fractional_Synthesis_Rate[i0, j] > fBest_Fractional_Synthesis_Rate_Time[i, j - 1] &&
							//fMinSlope > Fractional_Synthesis_Rate[i0, j] - fBest_Fractional_Synthesis_Rate_Time[i, j - 1])
						//if (fMinSlope > Fractional_Synthesis_Rate[i0, j] - fBest_Fractional_Synthesis_Rate_Time[i, j - 1] )
						{
							fMinSlope = Fractional_Synthesis_Rate[i0, j] - fBest_Fractional_Synthesis_Rate_Time[i, j - 1];

							Score_Fractional_Synthesis_Rate_Time[i, j] = fSum;

							index_X[i, j] = i0;

							index_Y[i, j] = j;

							fBest_Fractional_Synthesis_Rate_Time[i, j] = Fractional_Synthesis_Rate[i0, j];

						}



					}
				}
				//puts("");

			}
		}



	}

	for (i = 0; i < Nrow; i++)
	{
		for (j = 0; j < Ncol; j++)
		{
			printf("Score[%d, %d] = %5.2f ", i, j, Score_Fractional_Synthesis_Rate_Time[i, j]);
		}
		puts("");
	}

	printf("Best path\n");


	String^ temp_path = "";
	for (i = 0; i < Nrow; i++)
	{
		temp_path = "";
		for (j = 0; j < Ncol; j++)
		{
			printf("Next[%d, %d] = %d %d  ", i, j, index_X[i, j], index_Y[i, j]);
			temp_path = temp_path + i.ToString() + "," + j.ToString() + " ";
			//temp_path = temp_path + fBest_Fractional_Synthesis_Rate_Time[i, j].ToString() +" ";
		}
		puts("");
		bestpaths->Add(temp_path->Trim());
		bestpaths_and_score->Add(temp_path->Trim());
	}

	ftemp = -1000000.;

	for (i = 0; i < Nrow; i++)
	{
		printf("%10.5f\n", Score_Fractional_Synthesis_Rate_Time[i, Ncol - 1]);
		bestpaths_score->Add(Score_Fractional_Synthesis_Rate_Time[i, Ncol - 1]);
		bestpaths_and_score[i] = bestpaths_and_score[i] + "#" + (Score_Fractional_Synthesis_Rate_Time[i, Ncol - 1]).ToString();
		if (Score_Fractional_Synthesis_Rate_Time[i, Ncol - 1] > ftemp)
		{
			iMax = i;

			ftemp = Score_Fractional_Synthesis_Rate_Time[i, Ncol - 1];
		}
	}








	return bestpaths_and_score;
}


float Label_Path::Diff_Fraction_Rates(float FRate1, float FRate2)
{
	if (FRate2 >= FRate1)
	{
		if (false)
			if (FRate2 > 1.1)
			{
				return -1. * (FRate2 - 1.1);
			}
			else
			{
				return 1.0;
			}
		else
		{
			return 1.2 - (FRate2 - FRate1);
		}

	}
	else
	{
		if (FRate2 < FRate1)
		{
			return -1.0 * (FRate1 - FRate2);
		}
		//else
			//return -1.0;
	}
}




#include<stdio.h>
#include<math.h>
#include<stdlib.h>




int Kronecker_Delta(int i, int j)
{
if (i == j)
	return 1;
else
	return 0;
}


//**********************************************

int C(int j, int Number_of_Collocation_points)
{
return 1 + Kronecker_Delta(j,0) + Kronecker_Delta(j,Number_of_Collocation_points);
}

//**********************************************

double * Generate_Chebyshev_Collocation_Points(int number_of_points)
{

double * Collocation_Points = malloc(sizeof(double)*number_of_points);  

for(int c =0; c < number_of_points; c++)
{
	Collocation_Points[c] = cos(M_PI*c/number_of_points);
}
return Collocation_Points; 
}



int main()
{

double * t = Generate_Chebyshev_Collocation_Points(30);

printf("%lf \n",t[7]);

free(t);

}




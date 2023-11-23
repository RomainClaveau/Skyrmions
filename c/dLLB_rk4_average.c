#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <inttypes.h>
#include "ran3.c"

#define sqr(x) ((x)*(x))

#define pi 3.141592653589793238462643383279502884197169399375105820974
#define tau 1e-7
#define lambda 0.1
#define D 1.0
#define N 1

// Computing epsilon symbol e(i,j,k)
double epsilon(int i, int j, int k)
{
	int e[3] = {i, j, k};
	
	int cp1[3] = {1, 2, 3};
	int cp2[3] = {3, 1, 2};
	int cp3[3] = {2, 3, 1};
	
	if(i == cp1[0] && j == cp1[1] && k == cp1[2]) return 1.0;
	if(i == cp2[0] && j == cp2[1] && k == cp2[2]) return 1.0;
	if(i == cp3[0] && j == cp3[1] && k == cp3[2]) return 1.0;
	if(i == j || i == k || j == k) return 0.0;
	return -1.0;
}

// Function f(s_i) such that d<s_i>/dt = f(s_i)
double func_s(double t, int i, int j, int k, double * s, double * w, double S[3][3], double O[3][3], double add_time, double add_space)
{
	return (1.0 / (1.0 + sqr(lambda))) * (epsilon(i, j, k) * w[j] * (s[k] + add_space) + epsilon(i, j, k) * O[j][k] + lambda * w[i] * O[k][k] - lambda * w[k] * S[k][i]);
}

// Function f(O_ij) such that d<O_ij>/dt = f(O_ij)
double func_O(double t, int i, int j, int k, int l, double * s, double * w, double S[3][3], double O[3][3], double add_time, double add_space)
{
	return ((-1.0) * O[i][j] / tau) + (1.0 / (1.0 + sqr(lambda))) * (epsilon(j, k, l) * w[k] * (O[i][l] + add_space) + D * epsilon(j, i, l) * s[l] / tau  + 2.0 * lambda * w[j] * s[l] * (O[i][l] + add_space) - lambda * w[l] * (s[l] * O[i][j] + s[j] * (O[i][l] + add_space)));
}

// Function f(S_ij) such that d<S_ij>/dt = f(S_ij)
double func_S(double t, int i, int j, int k, int l, double * s, double * w, double S[3][3], double O[3][3], double add_time, double add_space)
{
	return (1.0 / (1.0 + sqr(lambda))) * (epsilon(j, k, l) * w[k] * (S[i][l] + add_space) + epsilon(j, k, l) * (s[i] * O[k][l] + s[l] * O[k][i]) + lambda * w[j] * (s[i] * (S[l][l] + add_space) + 2.0 * s[l] * (S[i][l] + add_space) - 2.0 * s[i] * sqr(s[l])) - lambda * w[l] * (s[i] * (S[l][j] + add_space) + s[l] * (S[i][j] + add_space) + s[j] * (S[i][l] + add_space) - 2.0 * s[i] * s[l] * s[j])) + (1.0 / (1.0 + sqr(lambda))) * (epsilon(i, k, l) * w[k] * (S[j][l] + add_space) + epsilon(i, k, l) * (s[j] * O[k][l] + s[l] * O[k][j]) + lambda * w[i] * (s[j] * (S[l][l] + add_space) + 2.0 * s[l] * (S[j][l] + add_space) - 2.0 * s[j] * sqr(s[l])) - lambda * w[l] * (s[j] * (S[l][i] + add_space) + s[l] * (S[j][i] + add_space) + s[i] * (S[j][l] + add_space) - 2.0 * s[j] * s[l] * s[i]));
}

// Main loop
int main(void)
{
	// Creating file
	remove("data");
	FILE * file = fopen("data", "a+");
	
	long seed = 123456789;
	
	double tmin = 0.0;
	double tmax = 1e-1;
	double dt = 1e-8;
	double t;
	int s_arr = 1e4;
	
	double average_s[s_arr][3];
	
	for(int i = 0; i < s_arr; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			average_s[i][j] = 0.0;
		}
	}
	
	for(int c = 0; c < N; c++)
	{
		t = tmin; 
			
		// Random values
		double r1 = ran3(&seed);
		double r2 = ran3(&seed);
		double r3 = ran3(&seed);
		double r4 = ran3(&seed);
		double r5 = ran3(&seed);
		double r6 = ran3(&seed);
		
		double rn1 = sqrt(-2.0 * log(r2)) * cos(2.0 * pi * r1);
		double rn2 = sqrt(-2.0 * log(r4)) * sin(2.0 * pi * r3);
		double rn3 = sqrt(-2.0 * log(r6)) * cos(2.0 * pi * r5);
		
		double s[3] = {1.0, 0.0, 0.0};
		double w[3] = {0.0, 0.0, pi/5.0};
		double o[3] = {rn1, rn2, rn3};
		double S[3][3]; // S_ij = <s_i s_j>
		double O[3][3]; // O_ij = <s_i o_j>
		
		double new_s[3] = {1.0, 0.0, 0.0};
		double new_O[3][3];
		double new_S[3][3];
		
		double k1;
		double k2;
		double k3;
		double k4;
		
		// Computing initial S_ij tensor
		for(int i = 0; i < 3; i++){for(int j = 0; j < 3; j++){S[i][j] = s[i] * s[j];}}
		
		// Computing initial O_ij tensor
		for(int i = 0; i < 3; i++){for(int j = 0; j < 3; j++){O[i][j] = o[i] * s[j];}}
		
		while(t < tmax)
		{
			for(int i = 0; i < 3; i++)
			{
				new_s[i] = 0.0;
				
				for(int j = 0; j < 3; j++)
				{
					new_O[i][j] = 0.0;
					new_S[i][j] = 0.0;
				}
			}
			
			for(int i = 0; i < 3; i++)
			{
				for(int j = 0; j < 3; j++)
				{
					for(int k = 0; k < 3; k++)
					{
						// Computing RK4 for <s_i>
						k1 = func_s(t, i, j, k, s, w, S, O, 0.0, 0.0);
						k2 = func_s(t, i, j, k, s, w, S, O, dt/2.0, dt * k1 / 2.0);
						k3 = func_s(t, i, j, k, s, w, S, O, dt/2.0, dt * k2 / 2.0);
						k3 = func_s(t, i, j, k, s, w, S, O, dt, dt * k3);
						
						new_s[i] += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
														
						for(int l = 0; l < 3; l++)
						{
							// Computing RK4 for <O_ij>
							k1 = func_O(t, i, j, k, l, s, w, S, O, 0.0, 0.0);
							k2 = func_O(t, i, j, k, l, s, w, S, O, dt/2.0, dt * k1 / 2.0);
							k3 = func_O(t, i, j, k, l, s, w, S, O, dt/2.0, dt * k2 / 2.0);
							k3 = func_O(t, i, j, k, l, s, w, S, O, dt, dt * k3);
							
							new_O[i][j] += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
							
							// Computing RK4 for <S_ij>
							k1 = func_S(t, i, j, k, l, s, w, S, O, 0.0, 0.0);
							k2 = func_S(t, i, j, k, l, s, w, S, O, dt/2.0, dt * k1 / 2.0);
							k3 = func_S(t, i, j, k, l, s, w, S, O, dt/2.0, dt * k2 / 2.0);
							k3 = func_S(t, i, j, k, l, s, w, S, O, dt, dt * k3);
							
							new_S[i][j] += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
						}
					}
				}
			}
			
			for(int i = 0; i < 3; i++)
			{
				if(s[i] > 1e5 || s[i] < -1e5) printf("error");
				s[i] += new_s[i];
				
				for(int j = 0; j < 3; j++)
				{
					O[i][j] += new_O[i][j];
					S[i][j] += new_S[i][j];
				}
			}
					
			int interval = (int)((t - tmin) / dt);
			
			if(interval % 1000 == 0)
			{				
				for(int i = 0; i < 3; i++)
				{
					average_s[interval / 1000][i] += s[i] / N;
				}
			}
			
			if(interval % 10000 == 0)
			{
				printf("%i\n", interval);
			}
				
			t += dt;
		}
	}
	
	for(int i = 0; i < s_arr; i++)
	{
		fprintf(file, "%i\t%g\t%g\t%g\n", i, average_s[i][0], average_s[i][1], average_s[i][2]);
	}
	
	return 0;
}

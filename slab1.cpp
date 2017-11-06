#include<iostream>
#include<fstream>
#include<math.h>
#include<stdio.h>
#include<cstdlib>
#include<cmath>

#define pi 3.14159265
#define NUMSIMULATIONS 5000000

using namespace std;

int main()
{

double R = 0;
double T = 0;
double A = 0;

double s;
double psi;

double ai, at;

double cost;
double sint;
double cosp;
double sinp;

double mu_p_x;
double mu_p_y;
double mu_p_z;

double w;
double del_w;

double r;

double del_r = 0.005;
double del_alpha = 1;
double del_omega;

double grid_r;
double grid_alpha;



double d = 0.02;

double g = 0.75;

double z0 = 0.0;
double z1 = d;

int ir;
int ialpha;

bool r_resolved;
bool alpha_resolved;

double mu_a = 15.0;
double mu_s = 90.0;
double mu_t = mu_a + mu_s;


int i;

srand(time(NULL));


double angle = 0;


for(i = 1; i < NUMSIMULATIONS; i = i + 1)
{
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	double mu_x = 0.0;
	double mu_y = 0.0;
	double mu_z = 1.0;

	w = 0.9;

	s = -log(((double)rand()/(RAND_MAX))) / mu_t;

	x = x + mu_x * s;
	y = y + mu_y * s;
	z = z + mu_z * s;


	bool escaped = false;

	while( !escaped )
	{
			s = -log(((double)rand()/(RAND_MAX))) / mu_t;	
			psi = 2.0 * pi * ((double)rand()/(RAND_MAX));
			cost = (0.5 / g) * (1.0 + pow(g, 2.0) - pow((1.0 - pow(g, 2.0)) / (1.0 - g + 2.0 * g * ((double)rand()/(RAND_MAX))), 2.0));
			/*cost = 2 * ((double)rand()/(RAND_MAX)) - 1;*/
			cosp = cos(psi);
			sint = sqrt(1.0 - pow(cost, 2.0));
		

			if (psi < pi)
			{
				sinp = sqrt(1.0 - pow(cosp, 2.0));
			}
			else
			{
				sinp = -sqrt(1.0 - pow(cosp, 2.0));
			}

			if (abs(mu_z) > 0.99999)
			{
				mu_p_x = sint * cosp;
				mu_p_y = sint * sinp;
			
				if ( mu_z >= 0)
				{
					mu_p_z = cost;
				}
				else
				{
					mu_p_z = -cost;
				}
			}
			else
			{
				mu_p_x = sint / sqrt(1.0 - pow(mu_z, 2.0)) * (mu_x * mu_z * cosp - mu_y * sinp) + mu_x * cost;  
                mu_p_y = sint / sqrt(1.0 - pow(mu_z, 2.0)) * (mu_y * mu_z * cosp + mu_x * sinp) + mu_y * cost;  
				mu_p_z = -sint * cosp * sqrt(1.0 - pow(mu_z, 2.0)) + mu_z * cost;	
			}

			if ( (mu_p_z < 0) && (z + mu_p_z * s < 0) )
			{
				s = (z0 - z) / mu_p_z;

				x = x + mu_p_x * s;
				y = y + mu_p_y * s;
				z = z + mu_p_z * s;
				
				r = sqrt( pow(x, 2.0) + pow(y, 2.0) );

				R = R + w;
				
				escaped = true;
			}
			else if ( (mu_p_z > 0) && ( (z + mu_p_z * s) > d) )
			{
				s = (z1 - z) / mu_p_z;

				x = x + mu_p_x * s;
				y = y + mu_p_y * s;
				z = z + mu_p_z * s;

				T = T + w;

				escaped = true;
			}
			else
			{
				x = x + mu_p_x * s;
				y = y + mu_p_y * s;
				z = z + mu_p_z * s;

				 
				double random = ((double)rand()/(RAND_MAX));
				if ( (w < 0.0001) && ( random < 0.1) )
				{
					del_w = w * mu_a / mu_t;
					w = w - del_w;
					A = A + del_w;
					w = 10 * w;
				}
				else if( (w < 0.0001) && ( random >= 0.1) ) 	
				{
					escaped = true;
					A = A + w;
				}
				else
				{
					del_w = w * mu_a / mu_t;
					w = w - del_w;
					A = A + del_w;
				}
		
				mu_x = mu_p_x;	
				mu_y = mu_p_y;
				mu_z = mu_p_z;
	
			}
	}
}





double Absorption = A/i;
double Reflectance = R/i;
double Transmittance = T/i;

double R_expected = 0.09734;
double Reflectance_error = ((R/i - R_expected)/(R_expected) * 100.0);

double T_expected = 0.66096;
double Transmittance_error = ((T/i - T_expected)/(T_expected) * 100.0);



cout << "Reflectance is "<< Reflectance << endl;
cout << "Transmittance is "<< Transmittance << endl;
cout << "Reflectance error is "<< Reflectance_error <<endl;
cout << "Transmittance error is "<< Transmittance_error << endl;



ofstream myfile ("R-a_angle.txt");
if ( myfile.is_open())
{
	for ( int k=1; k < 101; k++)
	{
		
		myfile << R_alpha [k] << "\r\n";
		
	}
	
}


	

return 0;

}











		

#include <iostream>
#include <stdio.h>

using namespace std;

typedef struct tagPART
{
	double x, y; //coordinates
	double vx, vy; //velocity
	int outFlag; //is particle inside the box (0 - yes, 1 - no)

} PART;


const int NC = 1000; //number of cells
const int NP = 10; //number of particles
const double PI = 3.14159265358979323846;
const double k = 1.38e-23; //Bolzman const, J/K
const double m = 1.6735575e-27; //mass of Hydrogen atom, kg
const double R = 2 * 120e-12; //Van-der-Vaals radius of hydrogen, m

double Temp = 300; //temperature, K
double L = NC; //box length in normalized units
double Vterm = sqrt(k*Temp/m); //termal speed, m/s
double tStep = 0.1; //time step, s

PART particles[NP];

double fmaxwell(double v, double Vt)			//Maxwell distribution-function
{
	return pow((m/(2*PI*k*Temp)), 3/2)*exp(-(m*v * v) / (Vt * Vt));
}

double GetRand(double xmin, double xmax)		//random distribution
{
	return (double)rand() / RAND_MAX * (xmax - xmin) + xmin;
}

double RandMaxwell(double Vterm)					//random value distributed by Maxwell
{
	double v = GetRand(-40, 40);
	double f = GetRand(0, 1);
	while (f > fmaxwell(v, Vterm))
	{
		v = GetRand(-40, 40);
		f = GetRand(0, 1);
	}
	return v;
}

class GasParticles {
	public:

		void MoveParticles(void)										//move all particles a step forward
		{
			double x, y, vx, vy, r;
			for (int i = 0; i < NP; i++)
			{
				x = particles[i].x + particles[i].vx * tStep;			//coordinates of particle
				y = particles[i].y + particles[i].vy * tStep;

				if ((x > L || x < 0 || y > L || y < 0) && particles[i].outFlag != 1)					//check if particle inside the box
				{
					particles[i].x = x;
					particles[i].y = y;

					FILE* fout = fopen("output.dat", "a");
					fprintf(fout, "%i %f %f %f %f\n", i, particles[i].x, particles[i].y, particles[i].vx, particles[i].vy);
					fclose(fout);
					particles[i].outFlag = 1;

				}
			}

			for (int i=0; i<NP; i++)
			{
				for (int j = 0; j < NP; j++)						//check for collision with another particle
				{
					r = sqrt(pow(particles[i].x - particles[j].x, 2) + pow(particles[i].y - particles[j].y, 2));	//distance between particles

					if (i != j && particles[i].outFlag != 1 && r <= R)	// if particle is not itself, inside the box and on sertain distance 
					{
						vx = -particles[i].vx;						//as a first approximation we assume that collisions are elastic
						vy = -particles[i].vy;

						particles[i].x = x;
						particles[i].y = y;
						particles[i].vx = vx;
						particles[i].vy = vy;
						particles[j].vx = -particles[j].vx;
						particles[j].vy = -particles[j].vy;

						break;										//only two-particle collision are considered

					}

					else
					{
						particles[i].x = x;
						particles[i].y = y;

					}
						
				}
					
				
			}
		}

		void Init(void)  //particle intialisation
		{
			int i;
			for (i = 0; i < NP; i++)
			{
				particles[i].x = GetRand(0, L);
				particles[i].y = 0.0;							//particle are appeared on a bottom of the box
				particles[i].vx = RandMaxwell(Vterm);
				particles[i].vy = RandMaxwell(Vterm);
				particles[i].outFlag = 0;
			}
		}
};

int main()
{
	srand(time(NULL));
	int T = 50;
	remove("output.dat");

	GasParticles Group1;
	Group1.Init();

	for (int t = 0; t < T; t++)
	{
		Group1.MoveParticles();
	}

	return 0;
}
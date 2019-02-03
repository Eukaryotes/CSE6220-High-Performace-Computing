#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <limits>
#include <iostream>
#include <iomanip>

#define PI 3.14159265
using namespace std;

int dboard(int N, int rank, int size){
	int times, inside=0, i;
	double r, theta;
	const double threshold = 1/sqrt(2);

	if(rank==0)
		times = N/size + N%size;
	else
		times = N/size;

	for(i=0; i<times; i++){
		r = 1.0f * rand()/RAND_MAX;
		theta = 360.0f * rand()/RAND_MAX;
		if(abs(sqrt(r)*cos(theta*PI/180))<=threshold && abs(sqrt(r)*sin(theta*PI/180))<=threshold)
			inside ++;
	}
	return inside;
}

int main(int argc, char **argv){
	
	MPI_Init(&argc, &argv);
	MPI_Comm comm = MPI_COMM_WORLD;

	int rank, size, N, R;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	int inside=0, M=0;
	double totalPI = 0.0f;

	if(rank==0){
		N = stoi(argv[1]);
		R = stoi(argv[2]);
	}

	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&R, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int round = R;
	while(round>0){
		inside = dboard(N, rank, size);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(&inside, &M, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank==0)
			totalPI += 2.0f*N/M;
		round--;
	}
	if(rank==0) cout << std::setprecision (15) << 1.0f * totalPI/R << endl;

	MPI_Finalize();
	return 0;
}

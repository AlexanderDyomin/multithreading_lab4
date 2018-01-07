#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include "mpi.h"
#include "ParallelSort.h"

using namespace std;

int main(int argc, char* argv[]) {
	if (argc != 2) {
		cout << "Name of input file needed";
		return -1;
	}

	MPI_Init(&argc, &argv);
	int procRank, procNum;
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	int *mas;
	int masSize;
	double start, finish;

	if (procRank == 0) {
		ifstream input(argv[1]);
		input >> masSize;
		mas = new int[masSize];
		for (int i = 0; i < masSize; ++i) {
			input >> mas[i];
		}
	}

	ofstream out("results.txt");
	double interval, avg;
	avg = 0;
	for (int i = 0; i < 10; ++i) {
		std::random_shuffle(mas, mas + masSize);
		start = MPI_Wtime();
		ParallelSort::sort(mas, masSize, 0);
		finish = MPI_Wtime();
		interval = finish - start;
		avg += interval;
		if (procRank == 0)
			out << interval << endl;
	}
	avg /= 10;
	if (procRank == 0)
		out << avg << endl;
	out.close();

	delete[] mas;
	MPI_Finalize();
	return 0;
}
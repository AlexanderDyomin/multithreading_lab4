#pragma once
class ParallelSort
{
public:
	ParallelSort();
	~ParallelSort();

	static void sort(int *mas, int masSize, int rootProc);

private:
	static void greaterThan(int *mas, int &masSize, int *masNew, int &masNewSize, int *masSend
		, int &masSendSize, int pivot);
	static void lessThan(int *mas, int &masSize, int *masNew, int &masNewSize, int *masSend
		, int &masSendSize, int pivot);
};


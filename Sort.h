#pragma once
#include <vector>

using namespace std;

namespace FW
{
typedef bool    (*SortCompareFunc) (void* data, int idxA, int idxB);    // Returns true if A should come before B.
typedef void    (*SortSwapFunc)    (void* data, int idxA, int idxB);    // Swaps A and B.

void sort(void* data, int start, int end, SortCompareFunc compareFunc, SortSwapFunc swapFunc, bool multicore = false);

template <class T> void sort(vector<T>& data, int start, int end, SortCompareFunc compareFunc, SortSwapFunc swapFunc, bool multicore = false);

//------------------------------------------------------------------------

template <class T> void sort(vector<T>& data, int start, int end, SortCompareFunc compareFunc, SortSwapFunc swapFunc, bool multicore)
{
    sort(&data[start], 0, end - start, compareFunc, swapFunc, multicore);
}

//------------------------------------------------------------------------
}
/*
class Sort
{
public:
	Sort(void);
	~Sort(void);
};
*/
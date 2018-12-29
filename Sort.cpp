#include "Sort.h"

using namespace FW;

//------------------------------------------------------------------------

#define QSORT_STACK_SIZE    32
#define QSORT_MIN_SIZE      16
#define MULTICORE_MIN_SIZE  (1 << 13)

//------------------------------------------------------------------------

namespace FW
{
static inline void  insertionSort   (int start, int size, void* data, SortCompareFunc compareFunc, SortSwapFunc swapFunc);
static inline int   median3         (int low, int high, void* data, SortCompareFunc compareFunc);
static int          partition       (int low, int high, void* data, SortCompareFunc compareFunc, SortSwapFunc swapFunc);
static void         qsort           (int low, int high, void* data, SortCompareFunc compareFunc, SortSwapFunc swapFunc);

}


//------------------------------------------------------------------------

void FW::insertionSort(int start, int size, void* data, SortCompareFunc compareFunc, SortSwapFunc swapFunc)
{
    //FW_ASSERT(compareFunc && swapFunc);
    //FW_ASSERT(size >= 0);

    for (int i = 1; i < size; i++)
    {
        int j = start + i - 1;
        while (j >= start && compareFunc(data, j + 1, j))
        {
            swapFunc(data, j, j + 1);
            j--;
        }
    }
}

//------------------------------------------------------------------------

int FW::median3(int low, int high, void* data, SortCompareFunc compareFunc)
{
    //FW_ASSERT(compareFunc);
    //FW_ASSERT(low >= 0 && high >= 2);

    int l = low;
    int c = (low + high) >> 1;
    int h = high - 2;

    if (compareFunc(data, h, l)) swap(l, h);
    if (compareFunc(data, c, l)) c = l;
    return (compareFunc(data, h, c)) ? h : c;
}

//------------------------------------------------------------------------

int FW::partition(int low, int high, void* data, SortCompareFunc compareFunc, SortSwapFunc swapFunc)
{
    // Select pivot using median-3, and hide it in the highest entry.

    swapFunc(data, median3(low, high, data, compareFunc), high - 1);

    // Partition data.

    int i = low - 1;
    int j = high - 1;
    for (;;)
    {
        do
            i++;
        while (compareFunc(data, i, high - 1));
        do
            j--;
        while (compareFunc(data, high - 1, j));

        //FW_ASSERT(i >= low && j >= low && i < high && j < high);
        if (i >= j)
            break;

        swapFunc(data, i, j);
    }

    // Restore pivot.

    swapFunc(data, i, high - 1);
    return i;
}

//------------------------------------------------------------------------

void FW::qsort(int low, int high, void* data, SortCompareFunc compareFunc, SortSwapFunc swapFunc)
{
    //FW_ASSERT(compareFunc && swapFunc);
    //FW_ASSERT(low <= high);

    int stack[QSORT_STACK_SIZE];
    int sp = 0;
    stack[sp++] = high;

    while (sp)
    {
        high = stack[--sp];
        //FW_ASSERT(low <= high);

        // Small enough or stack full => use insertion sort.

        if (high - low < QSORT_MIN_SIZE || sp + 2 > QSORT_STACK_SIZE)
        {
            insertionSort(low, high - low, data, compareFunc, swapFunc);
            low = high + 1;
            continue;
        }

        // Partition and sort sub-partitions.

        int i = partition(low, high, data, compareFunc, swapFunc);
        //FW_ASSERT(sp + 2 <= QSORT_STACK_SIZE);
        if (high - i > 2)
            stack[sp++] = high;
        if (i - low > 1)
            stack[sp++] = i;
        else
            low = i + 1;
    }
}

//------------------------------------------------------------------------

//------------------------------------------------------------------------

void FW::sort(void* data, int start, int end, SortCompareFunc compareFunc, SortSwapFunc swapFunc, bool multicore)
{
    //FW_ASSERT(start <= end);
    //FW_ASSERT(compareFunc && swapFunc);

    // Nothing to do => skip.

    if (end - start < 2)
        return;

    // Single-core.

    if (!multicore || end - start < MULTICORE_MIN_SIZE)
    {
    qsort(start, end, data, compareFunc, swapFunc);
}

    // Multicore.
	/*
    else
{
    TaskSpec* spec = new TaskSpec;
    spec->low = start;
    spec->high = end;
    spec->data = data;
    spec->compareFunc = compareFunc;
    spec->swapFunc = swapFunc;

    MulticoreLauncher launcher;
    MulticoreLauncher::Task task;
    task.launcher = &launcher;
    task.data = spec;
    qsortMulticore(task);
}*/
}

//------------------------------------------------------------------------


/*
Sort::Sort(void)
{
}


Sort::~Sort(void)
{
}
*/
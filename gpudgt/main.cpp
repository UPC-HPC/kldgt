#include "gpudgt.h"
#include "gpu_deghost3d.h"
#include "cpu_deghost3d.h"

int main (int argc, char **argv)
{
    cpu_deghost3d();
    gpu_deghost3d();
	return 0;
}


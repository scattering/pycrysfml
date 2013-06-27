#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct VectorN {
    int v[5];
} vectorN;

vectorN __struct_MOD_structret(void);

int main(int argc, void *argv[])
{
    vectorN retval = __struct_MOD_structret();
    printf("%d\n",retval.v[0]);
    return 0;
}


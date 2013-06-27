#include <string.h>
#include <stdlib.h>

typedef struct Vector {
    int x, y, z;
} vector;

typedef struct VectorN {
    int v[5];
} vectorN;


typedef struct {
    int Rot[3][3];
    float Tr[3];
} Sym_Oper_Type;

Sym_Oper_Type symop = { {0}, {0}};

Sym_Oper_Type make() {
    Sym_Oper_Type retval;
    symop.Rot[0][0] = 42;
    memcpy(&retval, &symop, sizeof(Sym_Oper_Type));
    symop.Rot[0][0] = 43;
    return retval;
}

Sym_Oper_Type* makep() {
    symop.Rot[0][0] = 67;
    return &symop;
}
Sym_Oper_Type* makepcopy() {
    Sym_Oper_Type *retval = malloc(sizeof(Sym_Oper_Type));
    symop.Rot[0][0] = 69;
    memcpy(retval, &symop, sizeof(Sym_Oper_Type));
    return retval;
}

void freecopy(Sym_Oper_Type *ptr)
{
    free(ptr);
}

int made(Sym_Oper_Type *ptr)
{
    return ptr->Rot[0][0];
}

vectorN structret()
{
    vectorN retval;
    retval.v[0] = 49;
    return retval;
}

vector add(vector* u, vector* v)
{
    vector retval;
    retval.x = u->x + v->x;
    retval.y = u->y + v->y;
    retval.z = u->z + v->z;
    return retval;
}

int triple(int a)
{
    return 3*a;
}

int norm(vector* v)
{
    return v->x*v->x + v->y*v->y + v->z*v->z;
}

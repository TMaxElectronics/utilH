#ifndef UTIL_INC
#define UTIL_INC

#define CONFIG_MAX_LINE_SIZE 128
#define CONFIG_MAX_LINE_COUNT 128

//a handy define to get the size of an array
#define arraySize(X) (sizeof(X) / sizeof(X[0]))

#if __has_include("ff.h")
#include "ff.h"

char * CONFIG_getKey(FIL * file, char * keyToFind);
#endif

int32_t PWL_getY(int32_t x, int32_t * pwl, uint32_t listSizeRows, uint32_t preComputedDerivative);
int32_t atoiFP(char * a, uint32_t strlen, int32_t baseExponent, uint32_t ignoreUnit);

uint32_t isAsciiNumber(char c);
uint32_t isAsciiSpecialCharacter(char c);

#endif
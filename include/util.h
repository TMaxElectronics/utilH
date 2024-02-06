#ifndef UTIL_INC
#define UTIL_INC

#define CONFIG_MAX_LINE_SIZE 128
#define CONFIG_MAX_LINE_COUNT 128

#include "ff.h"

char * CONFIG_getKey(FIL * file, char * keyToFind);
int32_t PWL_getY(int32_t x, int32_t * pwl, uint32_t listSizeRows, uint32_t preComputedDerivative);

uint32_t isAsciiSpecialCharacter(char c);

#endif
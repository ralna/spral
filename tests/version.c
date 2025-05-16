#include "spral.h"
#include <stdio.h>

int main(void) {

    char version[11];
    get_spral_version(version);
    printf("SPRAL version: %s\n", version);

    return 0;
}

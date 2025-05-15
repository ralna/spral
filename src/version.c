#include "spral.h"
#include <string.h>
#include <stdio.h>

void get_spral_version(char version[11]) {
    strcpy(version, SPRAL_VERSION);
}

int main(void) {

    char version[11];
    get_spral_version(version);
    printf("SPRAL version: %s\n", version);

    return 0;
}

#include "spral.h"
#include <stdio.h>

void get_spral_version(int *year, int *month, int *day) {
    *year = SPRAL_VERSION_YEAR;
    *month = SPRAL_VERSION_MONTH;
    *day = SPRAL_VERSION_DAY;
}

int main(void) {

    int year, month, day;
    get_spral_version(&year, &month, &day);
    printf("SPRAL version: %4d.%02d.%02d\n", year, month, day);

    return 0;
}

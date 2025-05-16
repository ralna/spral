#include "spral.h"
#include <string.h>

void get_spral_version(char version[11]) {
    strcpy(version, SPRAL_VERSION);
}

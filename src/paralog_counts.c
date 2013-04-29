#include <stdio.h>
#include <stdlib.h>
#include "iniparser.h"

int main(int argc, char*argv[]) {
    printf("Paralog count loader\n");

    if (argc != 2) {
        fprintf(stderr, "Usage: %s <configuration.ini>\n", argv[0]);
        return -1;
    }

    // Load the configuration file.
    dictionary* ini;
    ini = iniparser_load(argv[1]);
    if (ini == NULL) {
        fprintf(stderr, "Cannot open configuration file: %s\n", "mipcounter.ini");
        return -1;
    }

    const int total_paralogs = iniparser_getsecnkeys(ini, "paralog_counts");
    char** paralog_keys = iniparser_getseckeys(ini, "paralog_counts");
    printf("Found %i paralogs\n", total_paralogs);

    int paralog_counts[total_paralogs];
    int i;
    for (i = 0; i < total_paralogs; i++) {
        paralog_counts[i] = iniparser_getint(ini, paralog_keys[i], -1);
        printf("Paralog %i: %i\n", i, paralog_counts[i]);
    }

    iniparser_freedict(ini);

    return 0;
}

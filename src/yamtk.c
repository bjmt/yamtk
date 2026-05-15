/*
 *   yamtk: Yet Another Motif Toolkit
 *   Copyright (C) 2025  Benjamin Jean-Marie Tremblay
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "yamscan.h"
#include "yamshuf.h"
#include "yamdedup.h"
#include "yamenr.h"
#include "yamme.h"
#include "yamref.h"
#include "version.h"

static void help(void) {
    printf(
        "yamtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "Usage:  yamtk <subcommand> [options]\n"
        "Available subcommands:\n"
        "    shuf       DNA/RNA sequence shuffling\n"
        "    scan       DNA/RNA motif scanning\n"
        "    dedup      Deduplicate overlapping ranges\n"
        "    enr        Motif enrichment\n"
        "    me         De novo motif elicitation\n"
        "    ref        Refine a PWM against positive sequences\n"
        "    version    Print the version number and exit\n"
        "    help       Print this message and exit\n"
        "For subcommand usage, try: yamtk <subcommand> -h\n"
        , YAMTK_VERSION, YAMTK_YEAR
    );
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Error: Missing subcommand\n");
        help();
        return EXIT_FAILURE;
    }
    if (strcmp(argv[1], "version") == 0) {
        fprintf(stderr, "v%s\n", YAMTK_VERSION);
        return EXIT_SUCCESS;
    } else if (strcmp(argv[1], "help") == 0) {
        help();
        return EXIT_SUCCESS;
    } else if (strcmp(argv[1], "shuf") == 0) {
        return main_shuf(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "scan") == 0) {
        return main_scan(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "dedup") == 0) {
        return main_dedup(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "enr") == 0) {
        return main_enr(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "me") == 0) {
        return main_me(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "ref") == 0) {
        return main_ref(argc - 1, argv + 1);
    } else {
        fprintf(stderr, "Error: Unknown subcommand '%s'; try 'help' for usage\n", argv[1]);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


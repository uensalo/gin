/*
 * fmd: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * permutation_parser.h is part of fmd
 *
 * fmd is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * fmd is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef FMD_PERMUTATION_PARSER_H
#define FMD_PERMUTATION_PARSER_H
#include "fmd_vector.h"
#include <stdio.h>

static fmd_vector_t *permutation_parse(FILE* file) {
    fmd_vector_t *vec;
    fmd_vector_init(&vec, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    char line[256];
    while (fgets(line, sizeof(line), file)) {
        int_t val;
        val = strtoll(line, NULL, 10);
        fmd_vector_append(vec, (void*)val);
    }
    fmd_vector_fit(vec);
    return vec;
}

#endif //FMD_PERMUTATION_PARSER_H

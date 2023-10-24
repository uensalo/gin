/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * permutation_parser.h is part of gin
 *
 * gin is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * gin is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef GIN_PERMUTATION_PARSER_H
#define GIN_PERMUTATION_PARSER_H
#include "gin_vector.h"
#include <stdio.h>

static gin_vector_t *permutation_parse(FILE* file) {
    gin_vector_t *vec;
    gin_vector_init(&vec, GIN_VECTOR_INIT_SIZE, &prm_fstruct);
    char line[256];
    while (fgets(line, sizeof(line), file)) {
        int_t val;
        val = strtoll(line, NULL, 10);
        gin_vector_append(vec, (void*)val);
    }
    gin_vector_fit(vec);
    return vec;
}

#endif //GIN_PERMUTATION_PARSER_H

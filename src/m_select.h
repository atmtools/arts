#pragma once

#include <matpack.h>

void Select(Vector& needles,
            const Vector& haystack,
            const ArrayOfIndex& needleind);

void Select(AscendingGrid& needles,
            const AscendingGrid& haystack,
            const ArrayOfIndex& needleind);

void Select(Sparse& needles,
            const Sparse& haystack,
            const ArrayOfIndex& needleind);

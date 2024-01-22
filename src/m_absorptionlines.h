#pragma once

#include "absorptionlines.h"

void abs_linesMirroring(ArrayOfAbsorptionLines& abs_lines, const String& type);

void abs_linesPopulation(ArrayOfAbsorptionLines& abs_lines, const String& type);

void abs_linesNormalization(ArrayOfAbsorptionLines& abs_lines,
                            const String& type);

void abs_linesCutoff(ArrayOfAbsorptionLines& abs_lines,
                     const String& type,
                     const Numeric& x);

void abs_linesLinemixingLimit(ArrayOfAbsorptionLines& abs_lines,
                              const Numeric& x);

void abs_linesLineShapeType(ArrayOfAbsorptionLines& abs_lines,
                            const String& type);

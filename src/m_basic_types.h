#pragma once

#include <matpack.h>

void VectorInsertGridPoints(  // WS Generic Output:
    Vector& og,               // Output grid
    // WS Generic Input:
    const Vector& ingrid,  // Input grid
    const Vector& points);

void VectorGaussian(Vector& y,
                    const Vector& x,
                    const Numeric& x0,
                    const Numeric& si,
                    const Numeric& fwhm);

void MatrixGaussian(Matrix& Y,
                    const Vector& x_row,
                    const Numeric& x0_row,
                    const Numeric& si_row,
                    const Numeric& fwhm_row,
                    const Vector& x_col,
                    const Numeric& x0_col,
                    const Numeric& si_col,
                    const Numeric& fwhm_col);

void VectorLogSpace(Vector& x,
                    const Numeric& start,
                    const Numeric& stop,
                    const Numeric& step);

void VectorLinSpace(Vector& x,
                    const Numeric& start,
                    const Numeric& stop,
                    const Numeric& step);

void VectorNLinSpace(Vector& x,
                     const Index& n,
                     const Numeric& start,
                     const Numeric& stop);

void VectorNLogSpace(Vector& x,
                     const Index& n,
                     const Numeric& start,
                     const Numeric& stop);

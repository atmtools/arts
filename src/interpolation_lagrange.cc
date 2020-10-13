#include "interpolation_lagrange.h"

namespace Interpolation {
Vector interpweights(const Lagrange& dim0) {return dim0.lx;}

Vector dinterpweights(const Lagrange& dim0) {return dim0.dlx;}

Matrix interpweights(const Lagrange& dim0, const Lagrange& dim1) {
  Matrix out (dim0.lx.nelem(), dim1.lx.nelem());
  for (Index i=0; i<dim0.lx.nelem(); i++)
  for (Index j=0; j<dim1.lx.nelem(); j++)
    out(i, j) = dim0.lx[i] * dim1.lx[j];
  return out;
}

Matrix dinterpweights(const Lagrange& dim0, const Lagrange& dim1, Index dim) {
  Matrix out (dim0.lx.nelem(), dim1.lx.nelem());
  for (Index i=0; i<dim0.lx.nelem(); i++)
  for (Index j=0; j<dim1.lx.nelem(); j++)
    out(i, j) = (dim==0 ? dim0.dlx[i] : dim0.lx[i]) *
                (dim==1 ? dim1.dlx[j] : dim1.lx[j]);
  return out;
}

Tensor3 interpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2) {
  Tensor3 out (dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem());
  for (Index i=0; i<dim0.lx.nelem(); i++)
  for (Index j=0; j<dim1.lx.nelem(); j++)
  for (Index k=0; k<dim2.lx.nelem(); k++)
    out(i, j, k) = dim0.lx[i] * dim1.lx[j] * dim2.lx[k];
  return out;
}

Tensor3 dinterpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, Index dim) {
  Tensor3 out (dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem());
  for (Index i=0; i<dim0.lx.nelem(); i++)
  for (Index j=0; j<dim1.lx.nelem(); j++)
  for (Index k=0; k<dim2.lx.nelem(); k++)
    out(i, j, k) = (dim==0 ? dim0.dlx[i] : dim0.lx[i]) *
                   (dim==1 ? dim1.dlx[j] : dim1.lx[j]) *
                   (dim==2 ? dim2.dlx[k] : dim2.lx[k]);
  return out;
}
 
Tensor4 interpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3) {
  Tensor4 out (dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem(), dim3.lx.nelem());
  for (Index i=0; i<dim0.lx.nelem(); i++)
  for (Index j=0; j<dim1.lx.nelem(); j++)
  for (Index k=0; k<dim2.lx.nelem(); k++)
  for (Index l=0; l<dim3.lx.nelem(); l++)
    out(i, j, k, l) = dim0.lx[i] * dim1.lx[j] * dim2.lx[k] * dim3.lx[l];
  return out;
}

Tensor4 dinterpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, Index dim) {
  Tensor4 out (dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem(), dim3.lx.nelem());
  for (Index i=0; i<dim0.lx.nelem(); i++)
  for (Index j=0; j<dim1.lx.nelem(); j++)
  for (Index k=0; k<dim2.lx.nelem(); k++)
  for (Index l=0; l<dim3.lx.nelem(); l++)
    out(i, j, k, l) = (dim==0 ? dim0.dlx[i] : dim0.lx[i]) *
                      (dim==1 ? dim1.dlx[j] : dim1.lx[j]) *
                      (dim==2 ? dim2.dlx[k] : dim2.lx[k]) *
                      (dim==3 ? dim3.dlx[l] : dim3.lx[l]);
  return out;
}

Tensor5 interpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4) {
  Tensor5 out (dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem(), dim3.lx.nelem(), dim4.lx.nelem());
  for (Index i=0; i<dim0.lx.nelem(); i++)
  for (Index j=0; j<dim1.lx.nelem(); j++)
  for (Index k=0; k<dim2.lx.nelem(); k++)
  for (Index l=0; l<dim3.lx.nelem(); l++)
  for (Index m=0; m<dim4.lx.nelem(); m++)
    out(i, j, k, l, m) = dim0.lx[i] * dim1.lx[j] * dim2.lx[k] * dim3.lx[l] * dim4.lx[m];
  return out;
}

Tensor5 dinterpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, Index dim) {
  Tensor5 out (dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem(), dim3.lx.nelem(), dim4.lx.nelem());
  for (Index i=0; i<dim0.lx.nelem(); i++)
  for (Index j=0; j<dim1.lx.nelem(); j++)
  for (Index k=0; k<dim2.lx.nelem(); k++)
  for (Index l=0; l<dim3.lx.nelem(); l++)
  for (Index m=0; m<dim4.lx.nelem(); m++)
    out(i, j, k, l, m) = (dim==0 ? dim0.dlx[i] : dim0.lx[i]) *
                         (dim==1 ? dim1.dlx[j] : dim1.lx[j]) *
                         (dim==2 ? dim2.dlx[k] : dim2.lx[k]) *
                         (dim==3 ? dim3.dlx[l] : dim3.lx[l]) *
                         (dim==4 ? dim4.dlx[m] : dim4.lx[m]);
  return out;
}

Tensor6 interpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5) {
  Tensor6 out (dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem(), dim3.lx.nelem(), dim4.lx.nelem(), dim5.lx.nelem());
  for (Index i=0; i<dim0.lx.nelem(); i++)
  for (Index j=0; j<dim1.lx.nelem(); j++)
  for (Index k=0; k<dim2.lx.nelem(); k++)
  for (Index l=0; l<dim3.lx.nelem(); l++)
  for (Index m=0; m<dim4.lx.nelem(); m++)
  for (Index n=0; n<dim5.lx.nelem(); n++)
    out(i, j, k, l, m, n) = dim0.lx[i] * dim1.lx[j] * dim2.lx[k] * dim3.lx[l] * dim4.lx[m] * dim5.lx[n];
  return out;
}

Tensor6 dinterpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5, Index dim) {
  Tensor6 out (dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem(), dim3.lx.nelem(), dim4.lx.nelem(), dim5.lx.nelem());
  for (Index i=0; i<dim0.lx.nelem(); i++)
  for (Index j=0; j<dim1.lx.nelem(); j++)
  for (Index k=0; k<dim2.lx.nelem(); k++)
  for (Index l=0; l<dim3.lx.nelem(); l++)
  for (Index m=0; m<dim4.lx.nelem(); m++)
  for (Index n=0; n<dim5.lx.nelem(); n++)
    out(i, j, k, l, m, n) = (dim==0 ? dim0.dlx[i] : dim0.lx[i]) *
                            (dim==1 ? dim1.dlx[j] : dim1.lx[j]) *
                            (dim==2 ? dim2.dlx[k] : dim2.lx[k]) *
                            (dim==3 ? dim3.dlx[l] : dim3.lx[l]) *
                            (dim==4 ? dim4.dlx[m] : dim4.lx[m]) *
                            (dim==5 ? dim5.dlx[n] : dim5.lx[n]);
  return out;
}

Tensor7 interpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5, const Lagrange& dim6) {
  Tensor7 out (dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem(), dim3.lx.nelem(), dim4.lx.nelem(), dim5.lx.nelem(), dim6.lx.nelem());
  for (Index i=0; i<dim0.lx.nelem(); i++)
  for (Index j=0; j<dim1.lx.nelem(); j++)
  for (Index k=0; k<dim2.lx.nelem(); k++)
  for (Index l=0; l<dim3.lx.nelem(); l++)
  for (Index m=0; m<dim4.lx.nelem(); m++)
  for (Index n=0; n<dim5.lx.nelem(); n++)
  for (Index o=0; o<dim6.lx.nelem(); o++)
    out(i, j, k, l, m, n, o) = dim0.lx[i] * dim1.lx[j] * dim2.lx[k] * dim3.lx[l] * dim4.lx[m] * dim5.lx[n] * dim6.lx[o];
  return out;
}

Tensor7 dinterpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5, const Lagrange& dim6, Index dim) {
  Tensor7 out (dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem(), dim3.lx.nelem(), dim4.lx.nelem(), dim5.lx.nelem(), dim6.lx.nelem());
  for (Index i=0; i<dim0.lx.nelem(); i++)
  for (Index j=0; j<dim1.lx.nelem(); j++)
  for (Index k=0; k<dim2.lx.nelem(); k++)
  for (Index l=0; l<dim3.lx.nelem(); l++)
  for (Index m=0; m<dim4.lx.nelem(); m++)
  for (Index n=0; n<dim5.lx.nelem(); n++)
  for (Index o=0; o<dim6.lx.nelem(); o++)
    out(i, j, k, l, m, n, o) = (dim==0 ? dim0.dlx[i] : dim0.lx[i]) *
                               (dim==1 ? dim1.dlx[j] : dim1.lx[j]) *
                               (dim==2 ? dim2.dlx[k] : dim2.lx[k]) *
                               (dim==3 ? dim3.dlx[l] : dim3.lx[l]) *
                               (dim==4 ? dim4.dlx[m] : dim4.lx[m]) *
                               (dim==5 ? dim5.dlx[n] : dim5.lx[n]) *
                               (dim==6 ? dim6.dlx[o] : dim6.lx[o]);
  return out;
}

Numeric interp(const ConstVectorView yi, const ConstVectorView iw, const Lagrange& dim0) {
  const std::array<Index, 1> size{dim0.lx.nelem(), };
  std::array<Index, 1> ittr{0, };
  const auto start_ittr = ittr;
  Numeric out(0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    out += iw[ittr[0] + dim0.pos0] * yi[ittr[0]];
  }
  return out;
}

Numeric interp(const ConstMatrixView yi, const ConstMatrixView iw, const Lagrange& dim0, const Lagrange& dim1) {
  const std::array<Index, 2> size{dim0.lx.nelem(), dim1.lx.nelem(), };
  std::array<Index, 2> ittr{0, 0, };
  const auto start_ittr = ittr;
  Numeric out(0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      out += iw(ittr[0] + dim0.pos0, ittr[1] + dim1.pos0) * yi(ittr[0], ittr[1]);
    }
  }
  return out;
}

Numeric interp(const ConstTensor3View yi, const ConstTensor3View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2) {
  const std::array<Index, 3> size{dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem(), };
  std::array<Index, 3> ittr{0, 0, 0, };
  const auto start_ittr = ittr;
  Numeric out(0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        out += iw(ittr[0] + dim0.pos0, ittr[1] + dim1.pos0, ittr[2] + dim2.pos0) * yi(ittr[0], ittr[1], ittr[2]);
      }
    }
  }
  return out;
}

Numeric interp(const ConstTensor4View yi, const ConstTensor4View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3) {
  const std::array<Index, 4> size{dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem(), dim3.lx.nelem(), };
  std::array<Index, 4> ittr{0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Numeric out(0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          out += iw(ittr[0] + dim0.pos0, ittr[1] + dim1.pos0, ittr[2] + dim2.pos0, ittr[3] + dim3.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3]);
        }
      }
    }
  }
  return out;
}

Numeric interp(const ConstTensor5View yi, const ConstTensor5View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4) {
  const std::array<Index, 5> size{dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem(), dim3.lx.nelem(), dim4.lx.nelem(), };
  std::array<Index, 5> ittr{0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Numeric out(0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            out += iw(ittr[0] + dim0.pos0, ittr[1] + dim1.pos0, ittr[2] + dim2.pos0, ittr[3] + dim3.pos0, ittr[4] + dim4.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4]);
          }
        }
      }
    }
  }
  return out;
}

Numeric interp(const ConstTensor6View yi, const ConstTensor6View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5) {
  const std::array<Index, 6> size{dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem(), dim3.lx.nelem(), dim4.lx.nelem(), dim5.lx.nelem(), };
  std::array<Index, 6> ittr{0, 0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Numeric out(0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            for (ittr[5]=start_ittr[5]; ittr[5]<size[5]; ittr[5]++) {
              out += iw(ittr[0] + dim0.pos0, ittr[1] + dim1.pos0, ittr[2] + dim2.pos0, ittr[3] + dim3.pos0, ittr[4] + dim4.pos0, ittr[5] + dim5.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4], ittr[5]);
            }
          }
        }
      }
    }
  }
  return out;
}

Numeric interp(const ConstTensor7View yi, const ConstTensor7View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5, const Lagrange& dim6) {
  const std::array<Index, 7> size{dim0.lx.nelem(), dim1.lx.nelem(), dim2.lx.nelem(), dim3.lx.nelem(), dim4.lx.nelem(), dim5.lx.nelem(), dim6.lx.nelem(), };
  std::array<Index, 7> ittr{0, 0, 0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Numeric out(0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            for (ittr[5]=start_ittr[5]; ittr[5]<size[5]; ittr[5]++) {
              for (ittr[6]=start_ittr[6]; ittr[6]<size[6]; ittr[6]++) {
                out += iw(ittr[0] + dim0.pos0, ittr[1] + dim1.pos0, ittr[2] + dim2.pos0, ittr[3] + dim3.pos0, ittr[4] + dim4.pos0, ittr[5] + dim5.pos0, ittr[6] + dim6.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4], ittr[5], ittr[6]);
              }
            }
          }
        }
      }
    }
  }
  return out;
}

Vector interp(const ConstMatrixView yi, const ConstVectorView iw, const Lagrange& dim0, const std::array<Index, 1>& axis) {
  const std::array<Index, 2> size{axis[0] == 0 ? dim0.lx.nelem() : yi.nrows(), axis[0] == 1 ? dim0.lx.nelem() : yi.ncols(), };
  const std::array<Index, 2> count{0 + (axis[0] >= 0), 0 + (axis[0] >= 1), };
  std::array<Index, 2> ittr{0, 0, };
  const auto start_ittr = ittr;
  Vector out(size[0 + count[0]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      out[ittr[0 + count[0]]] += iw[ittr[axis[0]] + dim0.pos0] * yi(ittr[0], ittr[1]);
    }
  }
  return out;
}

Vector interp(const ConstTensor3View yi, const ConstMatrixView iw, const Lagrange& dim0, const Lagrange& dim1, const std::array<Index, 2>& axis) {
  const std::array<Index, 3> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : yi.npages(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : yi.nrows(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : yi.ncols(), };
  const std::array<Index, 3> count{0 + (axis[0] >= 0) + (axis[1] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2), };
  std::array<Index, 3> ittr{0, 0, 0, };
  const auto start_ittr = ittr;
  Vector out(size[0 + count[0]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        out[ittr[0 + count[0]]] += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0) * yi(ittr[0], ittr[1], ittr[2]);
      }
    }
  }
  return out;
}

Vector interp(const ConstTensor4View yi, const ConstTensor3View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const std::array<Index, 3>& axis) {
  const std::array<Index, 4> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : axis[2] == 0 ? dim2.lx.nelem() : yi.nbooks(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : axis[2] == 1 ? dim2.lx.nelem() : yi.npages(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : axis[2] == 2 ? dim2.lx.nelem() : yi.nrows(), axis[0] == 3 ? dim0.lx.nelem() : axis[1] == 3 ? dim1.lx.nelem() : axis[2] == 3 ? dim2.lx.nelem() : yi.ncols(), };
  const std::array<Index, 4> count{0 + (axis[0] >= 0) + (axis[1] >= 0) + (axis[2] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1) + (axis[2] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2) + (axis[2] >= 2), 0 + (axis[0] >= 3) + (axis[1] >= 3) + (axis[2] >= 3), };
  std::array<Index, 4> ittr{0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Vector out(size[0 + count[0]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          out[ittr[0 + count[0]]] += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0, ittr[axis[2]] + dim2.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3]);
        }
      }
    }
  }
  return out;
}

Vector interp(const ConstTensor5View yi, const ConstTensor4View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const std::array<Index, 4>& axis) {
  const std::array<Index, 5> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : axis[2] == 0 ? dim2.lx.nelem() : axis[3] == 0 ? dim3.lx.nelem() : yi.nshelves(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : axis[2] == 1 ? dim2.lx.nelem() : axis[3] == 1 ? dim3.lx.nelem() : yi.nbooks(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : axis[2] == 2 ? dim2.lx.nelem() : axis[3] == 2 ? dim3.lx.nelem() : yi.npages(), axis[0] == 3 ? dim0.lx.nelem() : axis[1] == 3 ? dim1.lx.nelem() : axis[2] == 3 ? dim2.lx.nelem() : axis[3] == 3 ? dim3.lx.nelem() : yi.nrows(), axis[0] == 4 ? dim0.lx.nelem() : axis[1] == 4 ? dim1.lx.nelem() : axis[2] == 4 ? dim2.lx.nelem() : axis[3] == 4 ? dim3.lx.nelem() : yi.ncols(), };
  const std::array<Index, 5> count{0 + (axis[0] >= 0) + (axis[1] >= 0) + (axis[2] >= 0) + (axis[3] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1) + (axis[2] >= 1) + (axis[3] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2) + (axis[2] >= 2) + (axis[3] >= 2), 0 + (axis[0] >= 3) + (axis[1] >= 3) + (axis[2] >= 3) + (axis[3] >= 3), 0 + (axis[0] >= 4) + (axis[1] >= 4) + (axis[2] >= 4) + (axis[3] >= 4), };
  std::array<Index, 5> ittr{0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Vector out(size[0 + count[0]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            out[ittr[0 + count[0]]] += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0, ittr[axis[2]] + dim2.pos0, ittr[axis[3]] + dim3.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4]);
          }
        }
      }
    }
  }
  return out;
}

Vector interp(const ConstTensor6View yi, const ConstTensor5View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const std::array<Index, 5>& axis) {
  const std::array<Index, 6> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : axis[2] == 0 ? dim2.lx.nelem() : axis[3] == 0 ? dim3.lx.nelem() : axis[4] == 0 ? dim4.lx.nelem() : yi.nvitrines(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : axis[2] == 1 ? dim2.lx.nelem() : axis[3] == 1 ? dim3.lx.nelem() : axis[4] == 1 ? dim4.lx.nelem() : yi.nshelves(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : axis[2] == 2 ? dim2.lx.nelem() : axis[3] == 2 ? dim3.lx.nelem() : axis[4] == 2 ? dim4.lx.nelem() : yi.nbooks(), axis[0] == 3 ? dim0.lx.nelem() : axis[1] == 3 ? dim1.lx.nelem() : axis[2] == 3 ? dim2.lx.nelem() : axis[3] == 3 ? dim3.lx.nelem() : axis[4] == 3 ? dim4.lx.nelem() : yi.npages(), axis[0] == 4 ? dim0.lx.nelem() : axis[1] == 4 ? dim1.lx.nelem() : axis[2] == 4 ? dim2.lx.nelem() : axis[3] == 4 ? dim3.lx.nelem() : axis[4] == 4 ? dim4.lx.nelem() : yi.nrows(), axis[0] == 5 ? dim0.lx.nelem() : axis[1] == 5 ? dim1.lx.nelem() : axis[2] == 5 ? dim2.lx.nelem() : axis[3] == 5 ? dim3.lx.nelem() : axis[4] == 5 ? dim4.lx.nelem() : yi.ncols(), };
  const std::array<Index, 6> count{0 + (axis[0] >= 0) + (axis[1] >= 0) + (axis[2] >= 0) + (axis[3] >= 0) + (axis[4] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1) + (axis[2] >= 1) + (axis[3] >= 1) + (axis[4] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2) + (axis[2] >= 2) + (axis[3] >= 2) + (axis[4] >= 2), 0 + (axis[0] >= 3) + (axis[1] >= 3) + (axis[2] >= 3) + (axis[3] >= 3) + (axis[4] >= 3), 0 + (axis[0] >= 4) + (axis[1] >= 4) + (axis[2] >= 4) + (axis[3] >= 4) + (axis[4] >= 4), 0 + (axis[0] >= 5) + (axis[1] >= 5) + (axis[2] >= 5) + (axis[3] >= 5) + (axis[4] >= 5), };
  std::array<Index, 6> ittr{0, 0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Vector out(size[0 + count[0]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            for (ittr[5]=start_ittr[5]; ittr[5]<size[5]; ittr[5]++) {
              out[ittr[0 + count[0]]] += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0, ittr[axis[2]] + dim2.pos0, ittr[axis[3]] + dim3.pos0, ittr[axis[4]] + dim4.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4], ittr[5]);
            }
          }
        }
      }
    }
  }
  return out;
}

Vector interp(const ConstTensor7View yi, const ConstTensor6View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5, const std::array<Index, 6>& axis) {
  const std::array<Index, 7> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : axis[2] == 0 ? dim2.lx.nelem() : axis[3] == 0 ? dim3.lx.nelem() : axis[4] == 0 ? dim4.lx.nelem() : axis[5] == 0 ? dim5.lx.nelem() : yi.nlibraries(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : axis[2] == 1 ? dim2.lx.nelem() : axis[3] == 1 ? dim3.lx.nelem() : axis[4] == 1 ? dim4.lx.nelem() : axis[5] == 1 ? dim5.lx.nelem() : yi.nvitrines(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : axis[2] == 2 ? dim2.lx.nelem() : axis[3] == 2 ? dim3.lx.nelem() : axis[4] == 2 ? dim4.lx.nelem() : axis[5] == 2 ? dim5.lx.nelem() : yi.nshelves(), axis[0] == 3 ? dim0.lx.nelem() : axis[1] == 3 ? dim1.lx.nelem() : axis[2] == 3 ? dim2.lx.nelem() : axis[3] == 3 ? dim3.lx.nelem() : axis[4] == 3 ? dim4.lx.nelem() : axis[5] == 3 ? dim5.lx.nelem() : yi.nbooks(), axis[0] == 4 ? dim0.lx.nelem() : axis[1] == 4 ? dim1.lx.nelem() : axis[2] == 4 ? dim2.lx.nelem() : axis[3] == 4 ? dim3.lx.nelem() : axis[4] == 4 ? dim4.lx.nelem() : axis[5] == 4 ? dim5.lx.nelem() : yi.npages(), axis[0] == 5 ? dim0.lx.nelem() : axis[1] == 5 ? dim1.lx.nelem() : axis[2] == 5 ? dim2.lx.nelem() : axis[3] == 5 ? dim3.lx.nelem() : axis[4] == 5 ? dim4.lx.nelem() : axis[5] == 5 ? dim5.lx.nelem() : yi.nrows(), axis[0] == 6 ? dim0.lx.nelem() : axis[1] == 6 ? dim1.lx.nelem() : axis[2] == 6 ? dim2.lx.nelem() : axis[3] == 6 ? dim3.lx.nelem() : axis[4] == 6 ? dim4.lx.nelem() : axis[5] == 6 ? dim5.lx.nelem() : yi.ncols(), };
  const std::array<Index, 7> count{0 + (axis[0] >= 0) + (axis[1] >= 0) + (axis[2] >= 0) + (axis[3] >= 0) + (axis[4] >= 0) + (axis[5] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1) + (axis[2] >= 1) + (axis[3] >= 1) + (axis[4] >= 1) + (axis[5] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2) + (axis[2] >= 2) + (axis[3] >= 2) + (axis[4] >= 2) + (axis[5] >= 2), 0 + (axis[0] >= 3) + (axis[1] >= 3) + (axis[2] >= 3) + (axis[3] >= 3) + (axis[4] >= 3) + (axis[5] >= 3), 0 + (axis[0] >= 4) + (axis[1] >= 4) + (axis[2] >= 4) + (axis[3] >= 4) + (axis[4] >= 4) + (axis[5] >= 4), 0 + (axis[0] >= 5) + (axis[1] >= 5) + (axis[2] >= 5) + (axis[3] >= 5) + (axis[4] >= 5) + (axis[5] >= 5), 0 + (axis[0] >= 6) + (axis[1] >= 6) + (axis[2] >= 6) + (axis[3] >= 6) + (axis[4] >= 6) + (axis[5] >= 6), };
  std::array<Index, 7> ittr{0, 0, 0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Vector out(size[0 + count[0]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            for (ittr[5]=start_ittr[5]; ittr[5]<size[5]; ittr[5]++) {
              for (ittr[6]=start_ittr[6]; ittr[6]<size[6]; ittr[6]++) {
                out[ittr[0 + count[0]]] += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0, ittr[axis[2]] + dim2.pos0, ittr[axis[3]] + dim3.pos0, ittr[axis[4]] + dim4.pos0, ittr[axis[5]] + dim5.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4], ittr[5], ittr[6]);
              }
            }
          }
        }
      }
    }
  }
  return out;
}

Matrix interp(const ConstTensor3View yi, const ConstVectorView iw, const Lagrange& dim0, const std::array<Index, 1>& axis) {
  const std::array<Index, 3> size{axis[0] == 0 ? dim0.lx.nelem() : yi.npages(), axis[0] == 1 ? dim0.lx.nelem() : yi.nrows(), axis[0] == 2 ? dim0.lx.nelem() : yi.ncols(), };
  const std::array<Index, 3> count{0 + (axis[0] >= 0), 0 + (axis[0] >= 1), 0 + (axis[0] >= 2), };
  std::array<Index, 3> ittr{0, 0, 0, };
  const auto start_ittr = ittr;
  Matrix out(size[0 + count[0]], size[1 + count[1]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        out(ittr[0 + count[0]], ittr[1 + count[1]]) += iw[ittr[axis[0]] + dim0.pos0] * yi(ittr[0], ittr[1], ittr[2]);
      }
    }
  }
  return out;
}

Matrix interp(const ConstTensor4View yi, const ConstMatrixView iw, const Lagrange& dim0, const Lagrange& dim1, const std::array<Index, 2>& axis) {
  const std::array<Index, 4> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : yi.nbooks(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : yi.npages(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : yi.nrows(), axis[0] == 3 ? dim0.lx.nelem() : axis[1] == 3 ? dim1.lx.nelem() : yi.ncols(), };
  const std::array<Index, 4> count{0 + (axis[0] >= 0) + (axis[1] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2), 0 + (axis[0] >= 3) + (axis[1] >= 3), };
  std::array<Index, 4> ittr{0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Matrix out(size[0 + count[0]], size[1 + count[1]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          out(ittr[0 + count[0]], ittr[1 + count[1]]) += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3]);
        }
      }
    }
  }
  return out;
}

Matrix interp(const ConstTensor5View yi, const ConstTensor3View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const std::array<Index, 3>& axis) {
  const std::array<Index, 5> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : axis[2] == 0 ? dim2.lx.nelem() : yi.nshelves(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : axis[2] == 1 ? dim2.lx.nelem() : yi.nbooks(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : axis[2] == 2 ? dim2.lx.nelem() : yi.npages(), axis[0] == 3 ? dim0.lx.nelem() : axis[1] == 3 ? dim1.lx.nelem() : axis[2] == 3 ? dim2.lx.nelem() : yi.nrows(), axis[0] == 4 ? dim0.lx.nelem() : axis[1] == 4 ? dim1.lx.nelem() : axis[2] == 4 ? dim2.lx.nelem() : yi.ncols(), };
  const std::array<Index, 5> count{0 + (axis[0] >= 0) + (axis[1] >= 0) + (axis[2] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1) + (axis[2] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2) + (axis[2] >= 2), 0 + (axis[0] >= 3) + (axis[1] >= 3) + (axis[2] >= 3), 0 + (axis[0] >= 4) + (axis[1] >= 4) + (axis[2] >= 4), };
  std::array<Index, 5> ittr{0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Matrix out(size[0 + count[0]], size[1 + count[1]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            out(ittr[0 + count[0]], ittr[1 + count[1]]) += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0, ittr[axis[2]] + dim2.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4]);
          }
        }
      }
    }
  }
  return out;
}

Matrix interp(const ConstTensor6View yi, const ConstTensor4View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const std::array<Index, 4>& axis) {
  const std::array<Index, 6> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : axis[2] == 0 ? dim2.lx.nelem() : axis[3] == 0 ? dim3.lx.nelem() : yi.nvitrines(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : axis[2] == 1 ? dim2.lx.nelem() : axis[3] == 1 ? dim3.lx.nelem() : yi.nshelves(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : axis[2] == 2 ? dim2.lx.nelem() : axis[3] == 2 ? dim3.lx.nelem() : yi.nbooks(), axis[0] == 3 ? dim0.lx.nelem() : axis[1] == 3 ? dim1.lx.nelem() : axis[2] == 3 ? dim2.lx.nelem() : axis[3] == 3 ? dim3.lx.nelem() : yi.npages(), axis[0] == 4 ? dim0.lx.nelem() : axis[1] == 4 ? dim1.lx.nelem() : axis[2] == 4 ? dim2.lx.nelem() : axis[3] == 4 ? dim3.lx.nelem() : yi.nrows(), axis[0] == 5 ? dim0.lx.nelem() : axis[1] == 5 ? dim1.lx.nelem() : axis[2] == 5 ? dim2.lx.nelem() : axis[3] == 5 ? dim3.lx.nelem() : yi.ncols(), };
  const std::array<Index, 6> count{0 + (axis[0] >= 0) + (axis[1] >= 0) + (axis[2] >= 0) + (axis[3] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1) + (axis[2] >= 1) + (axis[3] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2) + (axis[2] >= 2) + (axis[3] >= 2), 0 + (axis[0] >= 3) + (axis[1] >= 3) + (axis[2] >= 3) + (axis[3] >= 3), 0 + (axis[0] >= 4) + (axis[1] >= 4) + (axis[2] >= 4) + (axis[3] >= 4), 0 + (axis[0] >= 5) + (axis[1] >= 5) + (axis[2] >= 5) + (axis[3] >= 5), };
  std::array<Index, 6> ittr{0, 0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Matrix out(size[0 + count[0]], size[1 + count[1]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            for (ittr[5]=start_ittr[5]; ittr[5]<size[5]; ittr[5]++) {
              out(ittr[0 + count[0]], ittr[1 + count[1]]) += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0, ittr[axis[2]] + dim2.pos0, ittr[axis[3]] + dim3.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4], ittr[5]);
            }
          }
        }
      }
    }
  }
  return out;
}

Matrix interp(const ConstTensor7View yi, const ConstTensor5View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const std::array<Index, 5>& axis) {
  const std::array<Index, 7> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : axis[2] == 0 ? dim2.lx.nelem() : axis[3] == 0 ? dim3.lx.nelem() : axis[4] == 0 ? dim4.lx.nelem() : yi.nlibraries(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : axis[2] == 1 ? dim2.lx.nelem() : axis[3] == 1 ? dim3.lx.nelem() : axis[4] == 1 ? dim4.lx.nelem() : yi.nvitrines(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : axis[2] == 2 ? dim2.lx.nelem() : axis[3] == 2 ? dim3.lx.nelem() : axis[4] == 2 ? dim4.lx.nelem() : yi.nshelves(), axis[0] == 3 ? dim0.lx.nelem() : axis[1] == 3 ? dim1.lx.nelem() : axis[2] == 3 ? dim2.lx.nelem() : axis[3] == 3 ? dim3.lx.nelem() : axis[4] == 3 ? dim4.lx.nelem() : yi.nbooks(), axis[0] == 4 ? dim0.lx.nelem() : axis[1] == 4 ? dim1.lx.nelem() : axis[2] == 4 ? dim2.lx.nelem() : axis[3] == 4 ? dim3.lx.nelem() : axis[4] == 4 ? dim4.lx.nelem() : yi.npages(), axis[0] == 5 ? dim0.lx.nelem() : axis[1] == 5 ? dim1.lx.nelem() : axis[2] == 5 ? dim2.lx.nelem() : axis[3] == 5 ? dim3.lx.nelem() : axis[4] == 5 ? dim4.lx.nelem() : yi.nrows(), axis[0] == 6 ? dim0.lx.nelem() : axis[1] == 6 ? dim1.lx.nelem() : axis[2] == 6 ? dim2.lx.nelem() : axis[3] == 6 ? dim3.lx.nelem() : axis[4] == 6 ? dim4.lx.nelem() : yi.ncols(), };
  const std::array<Index, 7> count{0 + (axis[0] >= 0) + (axis[1] >= 0) + (axis[2] >= 0) + (axis[3] >= 0) + (axis[4] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1) + (axis[2] >= 1) + (axis[3] >= 1) + (axis[4] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2) + (axis[2] >= 2) + (axis[3] >= 2) + (axis[4] >= 2), 0 + (axis[0] >= 3) + (axis[1] >= 3) + (axis[2] >= 3) + (axis[3] >= 3) + (axis[4] >= 3), 0 + (axis[0] >= 4) + (axis[1] >= 4) + (axis[2] >= 4) + (axis[3] >= 4) + (axis[4] >= 4), 0 + (axis[0] >= 5) + (axis[1] >= 5) + (axis[2] >= 5) + (axis[3] >= 5) + (axis[4] >= 5), 0 + (axis[0] >= 6) + (axis[1] >= 6) + (axis[2] >= 6) + (axis[3] >= 6) + (axis[4] >= 6), };
  std::array<Index, 7> ittr{0, 0, 0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Matrix out(size[0 + count[0]], size[1 + count[1]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            for (ittr[5]=start_ittr[5]; ittr[5]<size[5]; ittr[5]++) {
              for (ittr[6]=start_ittr[6]; ittr[6]<size[6]; ittr[6]++) {
                out(ittr[0 + count[0]], ittr[1 + count[1]]) += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0, ittr[axis[2]] + dim2.pos0, ittr[axis[3]] + dim3.pos0, ittr[axis[4]] + dim4.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4], ittr[5], ittr[6]);
              }
            }
          }
        }
      }
    }
  }
  return out;
}

Tensor3 interp(const ConstTensor4View yi, const ConstVectorView iw, const Lagrange& dim0, const std::array<Index, 1>& axis) {
  const std::array<Index, 4> size{axis[0] == 0 ? dim0.lx.nelem() : yi.nbooks(), axis[0] == 1 ? dim0.lx.nelem() : yi.npages(), axis[0] == 2 ? dim0.lx.nelem() : yi.nrows(), axis[0] == 3 ? dim0.lx.nelem() : yi.ncols(), };
  const std::array<Index, 4> count{0 + (axis[0] >= 0), 0 + (axis[0] >= 1), 0 + (axis[0] >= 2), 0 + (axis[0] >= 3), };
  std::array<Index, 4> ittr{0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Tensor3 out(size[0 + count[0]], size[1 + count[1]], size[2 + count[2]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          out(ittr[0 + count[0]], ittr[1 + count[1]], ittr[2 + count[2]]) += iw[ittr[axis[0]] + dim0.pos0] * yi(ittr[0], ittr[1], ittr[2], ittr[3]);
        }
      }
    }
  }
  return out;
}

Tensor3 interp(const ConstTensor5View yi, const ConstMatrixView iw, const Lagrange& dim0, const Lagrange& dim1, const std::array<Index, 2>& axis) {
  const std::array<Index, 5> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : yi.nshelves(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : yi.nbooks(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : yi.npages(), axis[0] == 3 ? dim0.lx.nelem() : axis[1] == 3 ? dim1.lx.nelem() : yi.nrows(), axis[0] == 4 ? dim0.lx.nelem() : axis[1] == 4 ? dim1.lx.nelem() : yi.ncols(), };
  const std::array<Index, 5> count{0 + (axis[0] >= 0) + (axis[1] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2), 0 + (axis[0] >= 3) + (axis[1] >= 3), 0 + (axis[0] >= 4) + (axis[1] >= 4), };
  std::array<Index, 5> ittr{0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Tensor3 out(size[0 + count[0]], size[1 + count[1]], size[2 + count[2]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            out(ittr[0 + count[0]], ittr[1 + count[1]], ittr[2 + count[2]]) += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4]);
          }
        }
      }
    }
  }
  return out;
}

Tensor3 interp(const ConstTensor6View yi, const ConstTensor3View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const std::array<Index, 3>& axis) {
  const std::array<Index, 6> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : axis[2] == 0 ? dim2.lx.nelem() : yi.nvitrines(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : axis[2] == 1 ? dim2.lx.nelem() : yi.nshelves(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : axis[2] == 2 ? dim2.lx.nelem() : yi.nbooks(), axis[0] == 3 ? dim0.lx.nelem() : axis[1] == 3 ? dim1.lx.nelem() : axis[2] == 3 ? dim2.lx.nelem() : yi.npages(), axis[0] == 4 ? dim0.lx.nelem() : axis[1] == 4 ? dim1.lx.nelem() : axis[2] == 4 ? dim2.lx.nelem() : yi.nrows(), axis[0] == 5 ? dim0.lx.nelem() : axis[1] == 5 ? dim1.lx.nelem() : axis[2] == 5 ? dim2.lx.nelem() : yi.ncols(), };
  const std::array<Index, 6> count{0 + (axis[0] >= 0) + (axis[1] >= 0) + (axis[2] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1) + (axis[2] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2) + (axis[2] >= 2), 0 + (axis[0] >= 3) + (axis[1] >= 3) + (axis[2] >= 3), 0 + (axis[0] >= 4) + (axis[1] >= 4) + (axis[2] >= 4), 0 + (axis[0] >= 5) + (axis[1] >= 5) + (axis[2] >= 5), };
  std::array<Index, 6> ittr{0, 0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Tensor3 out(size[0 + count[0]], size[1 + count[1]], size[2 + count[2]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            for (ittr[5]=start_ittr[5]; ittr[5]<size[5]; ittr[5]++) {
              out(ittr[0 + count[0]], ittr[1 + count[1]], ittr[2 + count[2]]) += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0, ittr[axis[2]] + dim2.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4], ittr[5]);
            }
          }
        }
      }
    }
  }
  return out;
}

Tensor3 interp(const ConstTensor7View yi, const ConstTensor4View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const std::array<Index, 4>& axis) {
  const std::array<Index, 7> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : axis[2] == 0 ? dim2.lx.nelem() : axis[3] == 0 ? dim3.lx.nelem() : yi.nlibraries(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : axis[2] == 1 ? dim2.lx.nelem() : axis[3] == 1 ? dim3.lx.nelem() : yi.nvitrines(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : axis[2] == 2 ? dim2.lx.nelem() : axis[3] == 2 ? dim3.lx.nelem() : yi.nshelves(), axis[0] == 3 ? dim0.lx.nelem() : axis[1] == 3 ? dim1.lx.nelem() : axis[2] == 3 ? dim2.lx.nelem() : axis[3] == 3 ? dim3.lx.nelem() : yi.nbooks(), axis[0] == 4 ? dim0.lx.nelem() : axis[1] == 4 ? dim1.lx.nelem() : axis[2] == 4 ? dim2.lx.nelem() : axis[3] == 4 ? dim3.lx.nelem() : yi.npages(), axis[0] == 5 ? dim0.lx.nelem() : axis[1] == 5 ? dim1.lx.nelem() : axis[2] == 5 ? dim2.lx.nelem() : axis[3] == 5 ? dim3.lx.nelem() : yi.nrows(), axis[0] == 6 ? dim0.lx.nelem() : axis[1] == 6 ? dim1.lx.nelem() : axis[2] == 6 ? dim2.lx.nelem() : axis[3] == 6 ? dim3.lx.nelem() : yi.ncols(), };
  const std::array<Index, 7> count{0 + (axis[0] >= 0) + (axis[1] >= 0) + (axis[2] >= 0) + (axis[3] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1) + (axis[2] >= 1) + (axis[3] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2) + (axis[2] >= 2) + (axis[3] >= 2), 0 + (axis[0] >= 3) + (axis[1] >= 3) + (axis[2] >= 3) + (axis[3] >= 3), 0 + (axis[0] >= 4) + (axis[1] >= 4) + (axis[2] >= 4) + (axis[3] >= 4), 0 + (axis[0] >= 5) + (axis[1] >= 5) + (axis[2] >= 5) + (axis[3] >= 5), 0 + (axis[0] >= 6) + (axis[1] >= 6) + (axis[2] >= 6) + (axis[3] >= 6), };
  std::array<Index, 7> ittr{0, 0, 0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Tensor3 out(size[0 + count[0]], size[1 + count[1]], size[2 + count[2]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            for (ittr[5]=start_ittr[5]; ittr[5]<size[5]; ittr[5]++) {
              for (ittr[6]=start_ittr[6]; ittr[6]<size[6]; ittr[6]++) {
                out(ittr[0 + count[0]], ittr[1 + count[1]], ittr[2 + count[2]]) += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0, ittr[axis[2]] + dim2.pos0, ittr[axis[3]] + dim3.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4], ittr[5], ittr[6]);
              }
            }
          }
        }
      }
    }
  }
  return out;
}

Tensor4 interp(const ConstTensor5View yi, const ConstVectorView iw, const Lagrange& dim0, const std::array<Index, 1>& axis) {
  const std::array<Index, 5> size{axis[0] == 0 ? dim0.lx.nelem() : yi.nshelves(), axis[0] == 1 ? dim0.lx.nelem() : yi.nbooks(), axis[0] == 2 ? dim0.lx.nelem() : yi.npages(), axis[0] == 3 ? dim0.lx.nelem() : yi.nrows(), axis[0] == 4 ? dim0.lx.nelem() : yi.ncols(), };
  const std::array<Index, 5> count{0 + (axis[0] >= 0), 0 + (axis[0] >= 1), 0 + (axis[0] >= 2), 0 + (axis[0] >= 3), 0 + (axis[0] >= 4), };
  std::array<Index, 5> ittr{0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Tensor4 out(size[0 + count[0]], size[1 + count[1]], size[2 + count[2]], size[3 + count[3]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            out(ittr[0 + count[0]], ittr[1 + count[1]], ittr[2 + count[2]], ittr[3 + count[3]]) += iw[ittr[axis[0]] + dim0.pos0] * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4]);
          }
        }
      }
    }
  }
  return out;
}

Tensor4 interp(const ConstTensor6View yi, const ConstMatrixView iw, const Lagrange& dim0, const Lagrange& dim1, const std::array<Index, 2>& axis) {
  const std::array<Index, 6> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : yi.nvitrines(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : yi.nshelves(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : yi.nbooks(), axis[0] == 3 ? dim0.lx.nelem() : axis[1] == 3 ? dim1.lx.nelem() : yi.npages(), axis[0] == 4 ? dim0.lx.nelem() : axis[1] == 4 ? dim1.lx.nelem() : yi.nrows(), axis[0] == 5 ? dim0.lx.nelem() : axis[1] == 5 ? dim1.lx.nelem() : yi.ncols(), };
  const std::array<Index, 6> count{0 + (axis[0] >= 0) + (axis[1] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2), 0 + (axis[0] >= 3) + (axis[1] >= 3), 0 + (axis[0] >= 4) + (axis[1] >= 4), 0 + (axis[0] >= 5) + (axis[1] >= 5), };
  std::array<Index, 6> ittr{0, 0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Tensor4 out(size[0 + count[0]], size[1 + count[1]], size[2 + count[2]], size[3 + count[3]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            for (ittr[5]=start_ittr[5]; ittr[5]<size[5]; ittr[5]++) {
              out(ittr[0 + count[0]], ittr[1 + count[1]], ittr[2 + count[2]], ittr[3 + count[3]]) += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4], ittr[5]);
            }
          }
        }
      }
    }
  }
  return out;
}

Tensor4 interp(const ConstTensor7View yi, const ConstTensor3View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const std::array<Index, 3>& axis) {
  const std::array<Index, 7> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : axis[2] == 0 ? dim2.lx.nelem() : yi.nlibraries(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : axis[2] == 1 ? dim2.lx.nelem() : yi.nvitrines(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : axis[2] == 2 ? dim2.lx.nelem() : yi.nshelves(), axis[0] == 3 ? dim0.lx.nelem() : axis[1] == 3 ? dim1.lx.nelem() : axis[2] == 3 ? dim2.lx.nelem() : yi.nbooks(), axis[0] == 4 ? dim0.lx.nelem() : axis[1] == 4 ? dim1.lx.nelem() : axis[2] == 4 ? dim2.lx.nelem() : yi.npages(), axis[0] == 5 ? dim0.lx.nelem() : axis[1] == 5 ? dim1.lx.nelem() : axis[2] == 5 ? dim2.lx.nelem() : yi.nrows(), axis[0] == 6 ? dim0.lx.nelem() : axis[1] == 6 ? dim1.lx.nelem() : axis[2] == 6 ? dim2.lx.nelem() : yi.ncols(), };
  const std::array<Index, 7> count{0 + (axis[0] >= 0) + (axis[1] >= 0) + (axis[2] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1) + (axis[2] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2) + (axis[2] >= 2), 0 + (axis[0] >= 3) + (axis[1] >= 3) + (axis[2] >= 3), 0 + (axis[0] >= 4) + (axis[1] >= 4) + (axis[2] >= 4), 0 + (axis[0] >= 5) + (axis[1] >= 5) + (axis[2] >= 5), 0 + (axis[0] >= 6) + (axis[1] >= 6) + (axis[2] >= 6), };
  std::array<Index, 7> ittr{0, 0, 0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Tensor4 out(size[0 + count[0]], size[1 + count[1]], size[2 + count[2]], size[3 + count[3]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            for (ittr[5]=start_ittr[5]; ittr[5]<size[5]; ittr[5]++) {
              for (ittr[6]=start_ittr[6]; ittr[6]<size[6]; ittr[6]++) {
                out(ittr[0 + count[0]], ittr[1 + count[1]], ittr[2 + count[2]], ittr[3 + count[3]]) += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0, ittr[axis[2]] + dim2.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4], ittr[5], ittr[6]);
              }
            }
          }
        }
      }
    }
  }
  return out;
}

Tensor5 interp(const ConstTensor6View yi, const ConstVectorView iw, const Lagrange& dim0, const std::array<Index, 1>& axis) {
  const std::array<Index, 6> size{axis[0] == 0 ? dim0.lx.nelem() : yi.nvitrines(), axis[0] == 1 ? dim0.lx.nelem() : yi.nshelves(), axis[0] == 2 ? dim0.lx.nelem() : yi.nbooks(), axis[0] == 3 ? dim0.lx.nelem() : yi.npages(), axis[0] == 4 ? dim0.lx.nelem() : yi.nrows(), axis[0] == 5 ? dim0.lx.nelem() : yi.ncols(), };
  const std::array<Index, 6> count{0 + (axis[0] >= 0), 0 + (axis[0] >= 1), 0 + (axis[0] >= 2), 0 + (axis[0] >= 3), 0 + (axis[0] >= 4), 0 + (axis[0] >= 5), };
  std::array<Index, 6> ittr{0, 0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Tensor5 out(size[0 + count[0]], size[1 + count[1]], size[2 + count[2]], size[3 + count[3]], size[4 + count[4]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            for (ittr[5]=start_ittr[5]; ittr[5]<size[5]; ittr[5]++) {
              out(ittr[0 + count[0]], ittr[1 + count[1]], ittr[2 + count[2]], ittr[3 + count[3]], ittr[4 + count[4]]) += iw[ittr[axis[0]] + dim0.pos0] * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4], ittr[5]);
            }
          }
        }
      }
    }
  }
  return out;
}

Tensor5 interp(const ConstTensor7View yi, const ConstMatrixView iw, const Lagrange& dim0, const Lagrange& dim1, const std::array<Index, 2>& axis) {
  const std::array<Index, 7> size{axis[0] == 0 ? dim0.lx.nelem() : axis[1] == 0 ? dim1.lx.nelem() : yi.nlibraries(), axis[0] == 1 ? dim0.lx.nelem() : axis[1] == 1 ? dim1.lx.nelem() : yi.nvitrines(), axis[0] == 2 ? dim0.lx.nelem() : axis[1] == 2 ? dim1.lx.nelem() : yi.nshelves(), axis[0] == 3 ? dim0.lx.nelem() : axis[1] == 3 ? dim1.lx.nelem() : yi.nbooks(), axis[0] == 4 ? dim0.lx.nelem() : axis[1] == 4 ? dim1.lx.nelem() : yi.npages(), axis[0] == 5 ? dim0.lx.nelem() : axis[1] == 5 ? dim1.lx.nelem() : yi.nrows(), axis[0] == 6 ? dim0.lx.nelem() : axis[1] == 6 ? dim1.lx.nelem() : yi.ncols(), };
  const std::array<Index, 7> count{0 + (axis[0] >= 0) + (axis[1] >= 0), 0 + (axis[0] >= 1) + (axis[1] >= 1), 0 + (axis[0] >= 2) + (axis[1] >= 2), 0 + (axis[0] >= 3) + (axis[1] >= 3), 0 + (axis[0] >= 4) + (axis[1] >= 4), 0 + (axis[0] >= 5) + (axis[1] >= 5), 0 + (axis[0] >= 6) + (axis[1] >= 6), };
  std::array<Index, 7> ittr{0, 0, 0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Tensor5 out(size[0 + count[0]], size[1 + count[1]], size[2 + count[2]], size[3 + count[3]], size[4 + count[4]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            for (ittr[5]=start_ittr[5]; ittr[5]<size[5]; ittr[5]++) {
              for (ittr[6]=start_ittr[6]; ittr[6]<size[6]; ittr[6]++) {
                out(ittr[0 + count[0]], ittr[1 + count[1]], ittr[2 + count[2]], ittr[3 + count[3]], ittr[4 + count[4]]) += iw(ittr[axis[0]] + dim0.pos0, ittr[axis[1]] + dim1.pos0) * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4], ittr[5], ittr[6]);
              }
            }
          }
        }
      }
    }
  }
  return out;
}

Tensor6 interp(const ConstTensor7View yi, const ConstVectorView iw, const Lagrange& dim0, const std::array<Index, 1>& axis) {
  const std::array<Index, 7> size{axis[0] == 0 ? dim0.lx.nelem() : yi.nlibraries(), axis[0] == 1 ? dim0.lx.nelem() : yi.nvitrines(), axis[0] == 2 ? dim0.lx.nelem() : yi.nshelves(), axis[0] == 3 ? dim0.lx.nelem() : yi.nbooks(), axis[0] == 4 ? dim0.lx.nelem() : yi.npages(), axis[0] == 5 ? dim0.lx.nelem() : yi.nrows(), axis[0] == 6 ? dim0.lx.nelem() : yi.ncols(), };
  const std::array<Index, 7> count{0 + (axis[0] >= 0), 0 + (axis[0] >= 1), 0 + (axis[0] >= 2), 0 + (axis[0] >= 3), 0 + (axis[0] >= 4), 0 + (axis[0] >= 5), 0 + (axis[0] >= 6), };
  std::array<Index, 7> ittr{0, 0, 0, 0, 0, 0, 0, };
  const auto start_ittr = ittr;
  Tensor6 out(size[0 + count[0]], size[1 + count[1]], size[2 + count[2]], size[3 + count[3]], size[4 + count[4]], size[5 + count[5]], 0.0);
  for (ittr[0]=start_ittr[0]; ittr[0]<size[0]; ittr[0]++) {
    for (ittr[1]=start_ittr[1]; ittr[1]<size[1]; ittr[1]++) {
      for (ittr[2]=start_ittr[2]; ittr[2]<size[2]; ittr[2]++) {
        for (ittr[3]=start_ittr[3]; ittr[3]<size[3]; ittr[3]++) {
          for (ittr[4]=start_ittr[4]; ittr[4]<size[4]; ittr[4]++) {
            for (ittr[5]=start_ittr[5]; ittr[5]<size[5]; ittr[5]++) {
              for (ittr[6]=start_ittr[6]; ittr[6]<size[6]; ittr[6]++) {
                out(ittr[0 + count[0]], ittr[1 + count[1]], ittr[2 + count[2]], ittr[3 + count[3]], ittr[4 + count[4]], ittr[5 + count[5]]) += iw[ittr[axis[0]] + dim0.pos0] * yi(ittr[0], ittr[1], ittr[2], ittr[3], ittr[4], ittr[5], ittr[6]);
              }
            }
          }
        }
      }
    }
  }
  return out;
}
}

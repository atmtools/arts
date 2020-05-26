import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pyarts.catalogues import Sparse
import ctypes as c

class Block(object):
    """
    Block of a covariance matrix.

    A block in a covariance matrix is identified by its row and column
    block indices. The block indices define to which retrieval quantities
    the covariance matrix block corresponds to.


    Attributes:

        i(int): Row-block index of the block
        j(int): Column-block index of the block
        row_start(int): Row index of the left- and uppermost element of the
            block w.r.t to the full covariance matrix.
        column_start(int) Column index of the left- and uppermost element of
            the block w.r.t. to the full covariance matrix.
        inverse(bool): Flag that indicates whether this block is part of
            the normal part of the covariance matrix or its inverse.
        matrix(numpy.ndarray of scipy.sparse): The matrix that the block
            consists of.
    """

    @classmethod
    def __from_covariance_matrix_block_struct__(cls, s, inverse):
        """
        Create a block from a :class:`CovarianceMatrixBlockStruct`
        returned from the ARTS API.

        Paramters:
            s: The :class:`CovarianceMatrixBlockStruct` to create the
               block from.
            inverse: Flag that indicates whether the block belongs to
               the normal part of the covariance matrix or its inverse.

        Returns: The :class:`Block` object that represents the given
            :class:`CovarianceMatrixBlockStruct`
        """

        i, j   = list(s.indices)
        rs, cs = list(s.position)
        m, n   = list(s.dimensions)

        if not s.inner_ptr:
            matrix = np.ctypeslib.as_array(c.cast(s.ptr, c.POINTER(c.c_double)), (m, n))
        else:
            nnz = s.nnz
            data = np.ctypeslib.as_array(c.cast(s.ptr, c.POINTER(c.c_double)),
                                         (nnz,))
            row_indices = np.ctypeslib.as_array(s.inner_ptr, (nnz,))
            col_starts  = np.ctypeslib.as_array(s.outer_ptr, (m + 1,))
            matrix = sp.sparse.csr_matrix((data, row_indices, col_starts),
                                          shape = (m, n))
        return Block(i, j, rs, cs, inverse, matrix)

    def __init__(self, i, j, row_start, column_start, inverse, matrix):
        """
        Parameters:
            i(int): The row-block index of the covariance matrix block.
            j(int): The column-block index of the covariance matrix block.
            row_start(int): Row index of the left- and uppermost element in
                in the block.
            column_start(int): Column index of the left- and uppermost element
                in the block
            inverse(bool): Flag indicating whether the block is part of the
                inverse of the covariance matrix or not.
            matrix(np.ndarray or sp.sparse): The matrix of which the block
                consists.

        """
        self._i = i
        self._j = j
        self._row_start    = row_start
        self._column_start = column_start
        self._inverse = inverse
        self._matrix = matrix

    #
    # Read-only properties
    #

    @property
    def i(self):
        return self._i

    @property
    def j(self):
        return self._j

    @property
    def row_start(self):
        return self._row_start

    @property
    def column_start(self):
        return self._column_start

    @property
    def inverse(self):
        return self._inverse

    @property
    def matrix(self):
        return self._matrix

    def write_xml(self, xmlwriter, attr = None):
        """
        Serialize block and write to xml stream.

        Opens a new tag for the block and writes the matrix into
        it. Attributes of the block are saved as attributes of
        the newly create block.

        Parameters:
            xmlwriter: The xml stream to which to write the block.
            attr(dict): Additional attributes that should be added
                the tag that is created for the block.
        """
        if attr is None:
            attr = {}

        attr["row_index"] = self.i
        attr["column_index"] = self.j
        attr["row_start"] = self.row_start
        attr["column_start"] = self.column_start
        attr["row_extent"], attr["column_extent"] = self.matrix.shape
        attr["is_inverse"] = int(self.inverse)

        if sp.sparse.issparse(self.matrix):
            if not type(self.matrix) == Sparse:
                # why? because I can ...
                self.matrix.__class__.write_xml = Sparse.write_xml
            attr["type"] = "Sparse"
        else:
            attr["type"] = "Dense"

        xmlwriter.open_tag('Block', attr)
        xmlwriter.write_xml(self.matrix)
        xmlwriter.close_tag()

class CovarianceMatrix(object):
    """:class:`CovarianceMatrix` representing the ARTS group of the same name

    The CovarianceMatrix class is used to represent covariance matrices for
    OEM calculations in ARTS. It is represented as a block diagonal matrix
    where each block represents covariances between two retrieval quantities.

    Since covariance matrices are symmetric only block lying on or above
    the diagonal are actually stored. The covariance matrix class is
    designed to hold both, the covariance matrix and its inverse.
    """
    #
    # Class methods
    #

    @classmethod
    def __from_variable_value_struct__(cls, s):
        """
        Implements ARTS-API interface for returning objects from
        an ARTS workspace.


        """
        from pyarts.workspace.api import arts_api

        n_blocks = s.dimensions[0]
        n_inv_blocks = s.dimensions[1]

        blocks = []
        for i in range(n_blocks):
            bs = arts_api.get_covariance_matrix_block(s.ptr, i, False)
            b = Block.__from_covariance_matrix_block_struct__(bs, False)
            blocks += [b]

        inv_blocks = []
        for i in range(n_inv_blocks):
            bs = arts_api.get_covariance_matrix_block(s.ptr, i, True)
            b = Block.__from_covariance_matrix_block_struct__(bs, True)
            inv_blocks += [b]

        return CovarianceMatrix(blocks, inv_blocks)


    @classmethod
    def from_xml(cls, xmlelement):
        """Load a covariance matrix from an ARTS XML fiile.

        Returns:
           The loaded covariance matrix as :class:`CovarianceMatrix` object
        """
        n_blocks = xmlelement.get("n_blocks")

        blocks = []
        inv_blocks = []
        for b in list(xmlelement):

            i = b.get("row_index")
            j = b.get("column_index")
            row_start    = int(b.get("row_start"))
            column_start = int(b.get("column_start"))
            inverse = bool(int(b.get("is_inverse")))
            matrix = b[0].value()

            b = Block(i, j, row_start, column_start, inverse, matrix)
            if inverse:
                inv_blocks += [b]
            else:
                blocks += [b]
        return CovarianceMatrix(blocks, inv_blocks)

    def __init__(self, blocks, inverse_blocks = [], workspace = None):
        """
        Create a covariance matrix object.

        Parameters:
            blocks(list): List containing the blocks that make up the
                covariance matrix.
            inverse_blocks: Blocks that make up the inverse of the
                covariance. Can be provided to avoid computation of
                the inverse of the covariance matrix.
            workspace: :class:`Workspace` to associate the covariance
                matrix to.
        """
        self._blocks         = blocks
        self._inverse_blocks = inverse_blocks
        self._workspace      = None

    #
    # Read-only properties
    #

    @property
    def blocks(self):
        return self._blocks

    @property
    def inverse_blocks(self):
        return self._inverse_blocks

    @property
    def workspace(self):
        return self._workspace

    #
    # Serialization
    #

    def write_xml(self, xmlwriter, attr = None):
        """
        Implements typhon xml serialization interface.

        Parameters:
            xmlwriter: The xml stream to which to write the block.
            attr(dict): Additional attributes that should be added
                the tag that is created for the block.
        """

        if attr is None:
            attr = {}

        attr['n_blocks'] = len(self.blocks) + len(self.inverse_blocks)
        xmlwriter.open_tag('CovarianceMatrix', attr)

        for b in self.blocks:
            xmlwriter.write_xml(b)
        for b in self.inverse_blocks:
            xmlwriter.write_xml(b)

        xmlwriter.close_tag()

    def to_dense(self):
        """Conversion to dense representation.

        Converts the covariance matrix to a 2-dimensional numpy.ndarray.

        Returns:
            The covariance matrix as dense matrix.

        """
        m = max([b.row_start + b.matrix.shape[0] for b in self.blocks])
        n = max([b.column_start + b.matrix.shape[1] for b in self.blocks])
        mat = np.zeros((m, n))
        for b in self.blocks:
            m0 = b.row_start
            n0 = b.column_start
            dm = b.matrix.shape[0]
            dn = b.matrix.shape[1]
            if sp.sparse.issparse(b.matrix):
                mat[m0 : m0 + dm, n0 : n0 + dn] = b.matrix.toarray()
            else:
                mat[m0 : m0 + dm, n0 : n0 + dn] = b.matrix
        return mat

def plot_covariance_matrix(covariance_matrix, ax = None):
    """
    Plots a covariance matrix.

    Parameters:
        covariance_matrix(:class:`CovarianceMatrix`): The covariance matrix
            to plot
        ax(matplotlib.axes): An axes object into which to plot the
            covariance matrix.
    """

    if ax is None:
        ax = plt.gca()

    for b in covariance_matrix.blocks:
        y = np.arange(b.row_start, b.row_start + b.matrix.shape[0] + 1) - 0.5
        x = np.arange(b.column_start, b.column_start + b.matrix.shape[1] + 1) - 0.5
        ax.pcolormesh(x, y, np.array(b.matrix.toarray()))



    m = max([b.row_start + b.matrix.shape[0] for b in covariance_matrix.blocks])
    n = max([b.column_start + b.matrix.shape[1] for b in covariance_matrix.blocks])
    ax.set_xlim([-0.5, n + 0.5])
    ax.set_ylim([m + 0.5, -0.5])

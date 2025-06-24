Documentation
#############

Several parts of ARTS generater automatic documentation.

This document shows some solutions if you run into a problem.

This is an incomplete document.

Examples folder
===============

The documentation for the webpage is generated from the
structure and layout of the files
in the ``examples/`` directory.
This happens via the ``python/doc/gen_examples.py`` script,

There are several complicated rules involved
in the generation of the documentation.
The basic idea is that every ``rst``, ``py``,
and ``ipynb`` file in the ``examples/`` directory
is considered a documentation file.

.. note::
  The ``py`` and ``ipynb`` files are also considered tests and will run during the CI.

Here is a short list of rules:

1.  The name of the folder is used as the title of the documentation page.
2.  The sorting of all the content is alphanumerical.  Put ``1-*``, ``2-*``,
    etc. in front of the file names to control the order.

    a. Using ``N-`` leaves the numbering in the title and TOC.
       Using ``N.`` does not and removes the number.
       The ordering remains alphanumerical.

    b. If an ``rst`` files has the same stem name as a ``py`` file,
       the ``rst`` file is reordered before the ``py`` file.

    c. The file ``README.rst`` is always put first in the
       folder regardless of other names.

3.  The ``rst`` files are copied as are.
4.  The ``py`` files are put into sphinx code-blocks.
5.  The ``ipynb`` files are copied and run using a Sphinx extension.

    .. note::
      We only allow one ``ipynb`` file per folder.  This is because of our understand of how the Sphinx extension works.

6.  When multiple ``py`` files or ``rst`` files are present in a
    folder, their names will be used as subtitles in the documentation.

    a. The ``README.rst`` file does not get a subtitle.
       Please generate its subtitles in the file itself.

7.  Main title level in the generated documents is
    the ``====`` level and subtitles are ``----``.
    See the Sphinx documentation on header levels for more information.


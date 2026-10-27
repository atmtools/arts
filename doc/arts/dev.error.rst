Errors inside the C++ code
##########################

Reporting errors with `debug.h`
===============================


What does it achieve?
---------------------

Errors are reported by throwing exceptions.  The tools in `debug.h` give
consistent output that allows 1) easy debugging and 2) clear messages about
the problems to the user.

How is it implemented?
----------------------

These are the tools in `debug.h` that are important:

1. `arts::user_error`, a function that throws an error
2. `ARTS_METHOD_ERROR_CATCH`, a macro that adds the function name to errors
   passing through a function

The error messages are thrown as exceptions of the `std::runtime_error` type.

In debug builds (`Debug` and `RelWithDebInfo`), `arts::user_error` adds the
file name, function name, line number, and column number of the call to the
message.  In release builds, the message is shown as it is.

arts::user_error
----------------

This is used for errors caused by the user.  It never returns, so it may end
a function that would otherwise need a return value.

.. code-block:: c++

    arts::user_error("This is a test message, {} is not {}", 2, 1);

Checks are written as an `if` around the call:

.. code-block:: c++

    if (x.size() != y.size())
      arts::user_error("x has {} elements but y has {}", x.size(), y.size());

The arguments of the message are then only evaluated when the check fails, so
a check that passes costs no more than its condition.  Do not compute anything
for the message outside of the `if`.

A message without arguments is used as it is, so braces in it are printed
as they are written.

Each call with arguments formats its message with `std::format`.  Prefer
making these calls in source files rather than in header files, as they
increase compilation times.

ARTS_METHOD_ERROR_CATCH
-----------------------

This is used in a function try-catch block.

.. code-block:: c++

    void func(...) try {
    } ARTS_METHOD_ERROR_CATCH

An error that passes through the function is thrown again with the function
name added on the line above it, and the original error indented below it.
When functions that use it call each other, the result is a trace of the
functions that failed, followed by the error:

.. code-block:: text

    void outer(...)
      void inner(...)
        x has 3 elements but y has 4

Guide to a valid use of the formatting string
---------------------------------------------

The formatting string should be a valid format string for the `std::format`
function.  This means that you need to use `{}` for each of the arguments that
you want to insert.  You may add extra formatting options to the `{}`, these
depend on the type you are formatting.

The formatting string of `arts::user_error` must use every argument.  A call
with an argument that the string does not use does not compile.

.. code-block:: c++

    std::format("I have {} apple(s), {} orange(s), and {} argument(s)", 1, 2, 3);
    // Should give the string: "I have 1 apple(s), 2 orange(s), and 3 argument(s)"

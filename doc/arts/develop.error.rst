Errors inside the C++ code
##########################

Use of macros in `debug.h`
==========================

What does it achieve?
---------------------

We need to log and throw at errors.  This is done using the macros in `debug.h`.

These macros give consistent output that allows 1) easy debugging and 2) clear
messages about the problems to the user.

How is it implemented?
----------------------

The macros are defined in `debug.h` and are used throughout the ARTS code.

There are the macros that are important in there.

1. ARTS_ASSERT
2. ARTS_USER_ERROR
3. ARTS_USER_ERROR_IF
4. ARTS_METHOD_ERROR_CATCH

All of these use the `std::format` library to format the error message.

All of these adds relevant line numbers and file names to the error message.
If the compiler was modern enough, even the function signature is added.

ARTS_ASSERT
-----------

This is used for errors caused by the programmer.

It is only available when debug flags are set.
The error message will be thrown as an exception of the `std::logic_error` type.

Upon setting `-DARTS_ASSERT_USE_C`, the macro will use the C `assert` function.
This will cause the program to abort if the condition is not met.
The error message will not be printed in this mode.

Example usage:

.. code-block:: c++

    ARTS_ASSERT(1 == 2, "This is a test message, {} is not {}", 2, 1);

ARTS_USER_ERROR
---------------

This is used for unconditional errors caused by the user.

This gives an error message for user errors.  It is always available.
The error message will be thrown as an exception of the `std::runtime_error` type.

.. code-block:: c++

    ARTS_USER_ERROR("This is a test message, {} is not {}", 2, 1);

ARTS_USER_ERROR_IF
------------------

This is used for conditional errors caused by the user.

This gives a conditional error for user errors.  It is always available.
The error message will be thrown as an exception of the `std::runtime_error` type.

.. code-block:: c++

    ARTS_USER_ERROR(1 != 2 "This is a test message, {} is not {}", 2, 1);

ARTS_METHOD_ERROR_CATCH
-----------------------

This is used in a function try-catch block.  It gives the method name.

This gives an error message for method errors.  It is always available.

The error message will be thrown as an exception of the `std::runtime_error` type.

.. code-block:: c++

    void func(...) try {
    } ARTS_METHOD_ERROR_CATCH

Guide to a valid use of the formatting string
---------------------------------------------

The formatting string should be a valid format string for the `std::format` function.
This means that you need to use `{}` for each the arguments that you want to insert.
You may add extra formatting options to the `{}`, these depend on the type you are formatting.

.. code-block:: c++

    std::format("I have {} apple(s), {} orange(s), and {} argument(s)", 1, 2, 3);
    // Should give the string: "I have 1 apple(s), 2 orange(s), and 3 argument(s)"

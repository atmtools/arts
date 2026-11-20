Use of ``std::formatter<T>``
============================

What does it achieve?
---------------------

We use formatters to allow converting ARTS classes to strings consistently.
This has some main uses:

1. See the content of the types in python using ``__repr__``, ``__str__``, or ``print``.
   Using formatters instead of stream operators is not just faster, but also
   more consistent.  And it allows us to specialize the formatting for each
   object using the python ``__format__`` interface to control some aspects of
   how the output looks.
2. Calls to ``arts::format`` inside the C++ code can be used to format strings
   consistently.  This is especially useful for logging and debugging,
   and our experience using it in XML IO code shows that it is
   significantly faster than the old stream-based solution to
   text manipulation.

How is it implemented?
----------------------

Most ARTS classes should provide a ``std::formatter<T>`` specialization to allow
formatting of the class using ``arts::format``, ``std::format`` and their ilks.
The only exceptions are those groups that are marked as ``value_type`` in the
workspace groups definition.

This is done by
specialization of the ``std::formatter<T>`` class template for the class in
question. The following is an outline of meta-code
for the specialization of this class that
must exist for most ARTS classes:

.. code-block:: cpp
  :caption: Meta code for formatting ARTS classes to strings
  :linenos:

  template <>
  struct std::formatter<ARTSTYPE> {
    MetaData...;

    [[nodiscard]] constexpr auto& inner_fmt();

    [[nodiscard]] constexpr const auto& inner_fmt() const;

    constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx);

    template <class FmtContext> FmtContext::iterator format(const ARTSTYPE& v, FmtContext& ctx) const;
  };

This specialization must live on the top-level namespace (otherwise it does not work).

The only two functions that are actually needed for the class to be
formattable are ``format`` and ``parse``.  The rest are helper functions and data that
allow us to make the format across multiple ARTS classes compatible.
The implementation of ``format`` is also specific to the class in question,
so its potential implementations will not be discussed here in detail,
but examples and best-practices will be given.

The ``MetaData`` is a placeholder for additional data that the class needs to store
to be formatted.  For sake of sanity in the code, it is recommended that this data
is either another ``std::formatter<T> fmt;`` or a ``format_tags tags;`` object, as explicitly named.
Any other names or more options are discouraged as it will break cross-class compatibility.
Only these two options are discussed here.

If your class's ``MetaData`` is ``std::formatter<T> fmt;`` then the following is the recommended implementation:

.. code-block:: cpp
  :caption: Meta code for formatting ARTS classes that inherit a formatter
  :linenos:

  template <>
  struct std::formatter<ARTSTYPE> {
    std::formatter<T> fmt;

    [[nodiscard]] constexpr auto& inner_fmt() {return fmt.inner_fmt();}

    [[nodiscard]] constexpr const auto& inner_fmt() const {return fmt.inner_fmt();}

    constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {return fmt.inner_fmt().parse(ctx);}

    template <class FmtContext> FmtContext::iterator format(const ARTSTYPE& v, FmtContext& ctx) const;
  };

If your class's ``MetaData`` is ``format_tags tags;`` then the following is the recommended implementation:

.. code-block:: cpp
  :caption: Meta code for formatting ARTS classes that owns a formatter
  :linenos:

  template <>
  struct std::formatter<ARTSTYPE> {
    format_tags tags;

    [[nodiscard]] constexpr auto& inner_fmt() {return *this;}

    [[nodiscard]] constexpr const auto& inner_fmt() const {return *this;}

    constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) { return parse_format_tags(tags, ctx); }

    template <class FmtContext> FmtContext::iterator format(const ARTSTYPE& v, FmtContext& ctx) const;
  };

The above ensures that there exists at most one ``tags`` object for each class, and that
the ``tags`` object is the only object that needs to be passed around to ensure compatibility.

Formatting in the C++ code
--------------------------

Use ``arts::format`` instead of ``std::format`` to create strings in the C++
code, and ``arts::print`` and ``arts::println`` instead of ``std::print`` and
``std::println`` to write them to ``stdout``, a ``std::FILE*`` or a
``std::ostream``.  They take the same arguments and check the format string at
compile time in the same way:

.. code-block:: cpp
  :caption: Formatting in the C++ code

  const std::string s = arts::format("{} has {} elements", name, v.size());
  arts::println(os, "{:IO}", v);

These format with ``arts::vformat`` and ``arts::vprint``, and
``arts::vformat_to`` writes to a ``std::format_context``.  All are compiled
once, in ``format_tags.cpp``.  ``std::format``, ``std::vformat``,
``std::format_to``, ``std::print`` and ``std::println`` instead compile the
whole formatting machinery of the standard library in every file that calls
them, which makes the build considerably slower.

Write formatted text to a stream with ``arts::print(os, ...)`` or
``arts::println(os, ...)``, not with ``os << arts::format(...)``.

Inside a formatter, ``tags.format`` writes builtin values, such as ``Numeric``
and ``Index``, with ``format_tags::single_format``.  These have specializations
for ``std::format_context`` in ``format_tags.h`` that are also compiled once.
A formatter that writes another builtin type, such as ``float``, should add a
specialization for it next to the others.

What formatter options are available?
-------------------------------------

The following options are available for the ``format_tags`` object:

1. ``bracket``. Activated by the ``"B"`` character in the format string.
   Will add brackets around the output if applicable.
2. ``short_str``. Activated by the ``"s"`` character in the format string.
   Will use a short string representation if applicable.
3. ``comma``. Activated by the ``","`` character in the format string.
   Will comma-separate elements if applicable.
4. ``names``. Activated by the ``"N"`` character in the format string.
   Will show the names of the elements if applicable.
5. ``io``. Activated by the ``"IO"`` characters in the format string.
   Will use the IO formatting so that XML streams work as expected.
6. ``q``. Activated by the ``"q"`` character in the format string.
   Will quote the quotable output.
7. ``n``. Activated by the ``"n"`` character in the format string.
   Switch the default separator behavior from space to newline (or vice versa).
8. ``depth``.  Cannot be activated, should be set manually in the
   ``parse`` method.

The default formatting string given to ``__str__`` is ``"{:qNB,}"``
and the default formatting string given to ``__repr__`` is ``"{:nsqNB,}"``.

What a type will do with these options is up to the type itself.  Generally,
a ``bracket`` type should add brackets around the output, a ``short_str`` type should
show a short string representation, a ``comma`` type should separate elements with commas,
a ``names`` type should show the names of the elements, an ``io`` type should use the IO
formatting, and a ``q`` type should quote the quotable output.  If you need the formatted
output to look a certain way, you need to write a test that enforces this.  The way these
options interact is also not defined, so combining them may lead to unexpected results.

The special ``depth`` option is a ``Size`` that is used to indicate how "deep" the
variable should be indented in the output.  It is used for instance in the ``matpack``
formatter to indicate how many levels of indentation should be used when formatting
the objects.

What do you need to think about when implementing a formatter?
--------------------------------------------------------------

1. Calling the ``format_tags`` object's ``format`` method is an efficient way
   to chain formatting calls together

  .. code-block:: cpp
    :caption: Example of chaining formatting calls
    :linenos:

    template <class FmtContext>
    FmtContext::iterator format(const ARTSTYPE& v, FmtContext& ctx) const {
      const std::string_view sep = tags.sep();
      return tags.format(ctx, v.m1, sep, v.m2);
    }

2. Whenever you format in a ``const char *``, that is anything in C++ that is directly
   written ``"I am a const char *"``, there will be a resulting ``'\0'`` character included
   in the formatted string.  This will cause problems if you intend to copy-paste the
   screen output, as the ``'\0'`` character will not be visible but there anyway.
   Sometimes the output will even be truncated at the first ``'\0'`` character.
   To avoid this, use ``std::string_view`` instead of ``const char *`` whenever possible.
   The easiest way to this is simply to write ``"I am a string_view"sv``, as the ``sv``
   makes it a ``std::string_view`` and avoids copying the last ``'\0'`` character.

Use of std::formatter<T>
========================

What does it achieve?
---------------------

We use formatters to allow converting ARTS classes to strings consistently.
This has some main uses:

1. See the content of the types in python using `__repr__`, `__str__`, or `print`.
   Using formatters instead of stream operators is not just faster, but also
   more consistent.  And it allows us to specialize the formatting for each
   object using the python `__format__` interface to control some aspects of
   how the output looks.
2. Calls to `std::format` inside the C++ code can be used to format strings
   consistently.  This is especially useful for logging and debugging,
   and may even be used in the XML IO code when `std::print` is available.

How is it implemented?
----------------------

Most ARTS classes should provide a `std::formatter<T>` specialization to allow
formatting of the class using `fmt::format` and its ilks.
The only exceptions are those groups that are marked as `value_type` in the
workspace groups definition.

This is done by
specialization of the `std::formatter<T>` class template for the class in
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
formattable are `format` and `parse`.  The rest are helper functions and data that
allow us to make the format across multiple ARTS classes compatible.
The implementation of `format` is also specific to the class in question,
so its potential implementations will not be discussed here in details,
but examples and best-practices will be given.

The `MetaData` is a placeholder for additional data that the class needs to store
to be formatted.  For sake of sanity in the code, it is recommended that this data
is either another `std::formatter<T> fmt;` or a `format_tags tags;` object, as explicitly named.
Any other names or more options are discouraged as it will break cross-class compatibility.
Only these two options are discussed here.

If your class's `MetaData` is `std::formatter<T> fmt;` then the following is the recommended implementation:

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

If your class's `MetaData` is `format_tags tags;` then the following is the recommended implementation:

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

The above ensures that there exists at most one `tags` object for each class, and that
the `tags` object is the only object that needs to be passed around to ensure compatibility.

What formatter options are available?
-------------------------------------

The following options are available for the `format_tags` object:

1. `bracket`. Activated by the `B` character in the format string.
2. `short_str`. Activated by the `s` character in the format string.
3. `comma`. Activated by the `,` character in the format string.
4. `names`. Activated by the `N` character in the format string.

The default formatting string given to `__str__` is `{:NB,}` and the default
formatting string given to `__repr__` is `{:sNB,}`.

What a type will do with these options is up to the type itself.

What do you need to think about when implementing a formatter?
--------------------------------------------------------------

1. Calling the `format_tags` object's `format` method is an efficient way
   to chain formatting calls together

  .. code-block:: cpp
    :caption: Example of chaining formatting calls
    :linenos:

    template <class FmtContext>
    FmtContext::iterator format(const ARTSTYPE& v, FmtContext& ctx) const {
      const std::string_view sep = tags.sep();
      return tags.format(ctx, v.m1, sep, v.m2);
    }

2. Whenever you format in a `const char *`, that is anything in C++ that is directly
   written `"I am a const char *"`, there will be a resulting `'\0'` character included
   in the formatted string.  This will cause problems if you intend to copy-paste the 
   screen output, as the `'\0'` character will not be visible but there anyway.
   To avoid this, use `std::string_view` instead of `const char *` whenever possible.
   The easiest way to this is simply to write `"I am a string_view"sv`, as the `sv`
   makes it a `std::string_view` and avoids copying the last `'\0'` character.

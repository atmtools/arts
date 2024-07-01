Use of std::formatter<T>
========================

Most ARTS classes should provide a `std::formatter<T>` specialization to allow
formatting of the class using `fmt::format` and its ilks.
The only exceptions are `String`, `Numeric`, and `Index` as we want to rely
on the most standard method of printing these values.

If you are

 This is done by
specialization of the `std::formatter<T>` class template for the class in
question. The following is an outline of meta-code
for the specialization of this class that
must exist for most ARTS classes:

.. code-block:: cpp
  :caption: src/core/matpack/rational.h
  :linenos:

  template <>
  struct std::formatter<ARTSTYPE> {
    format_tags tags;
    OR
    std::formatter<OTHER_ARTSTYPE> fmt;

    [[nodiscard]] constexpr auto& inner_fmt();
    
    template <typename T> void compat(std::formatter<T>& x) const;
    
    template <typename... Ts> constexpr void make_compat(std::formatter<Ts>&... xs);
    
    constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx);

    template <class FmtContext> FmtContext::iterator format(const ARTSTYPE& v, FmtContext& ctx) const;
  };

The only two functions that are actually needed for the class to be
formattable are `format` and `parse`.  The other functions are ARTS-specific.

The `METADATA...` must be defined class-by-class and made compatible with the other types.
There are two main ways to specialize a `std::formatter<T>` for your new ARTS class.  The first is to
have its metadata be a pure `format_tags tags;` data object.  One class that does this is `Rational`.
The second is to have the metadata be a single other `std::formatter<T>`.  An example for this 
is `Vector` and all owning `matpack_data` objects.

The `inner_fmt` function should return a reference to the `std::formatter<T>` that is responsible
for the formatting of the class.  If the `METADATA` is a `format_tags` object, then this function
returns `*this`, and if `METADATA` is another `std::formatter<T>`, then this function returns
that `std::formatter<T>` object.

The `compat` and `make_compat` functions are used to make the `std::formatter<T>` compatible with other `std::formatter<T>`s.
`compat` makes the 

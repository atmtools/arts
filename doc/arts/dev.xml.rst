XML IO
======

This is the main way that we perform file IO in ARTS.

Each workspace group must support XML file IO.
If a group is added that does not support XML file IO, ARTS cannot compile.
Support is done by overloading a templated struct, ``xml_io_stream``,
which defines the core concepts involved in the workflow.

How the group handles XML file IO inside this is up to the group itself.
Generally, all lower level overloads should already be defined, so it
is common to just chain several calls to other XML IO methods.

The ``xml_io_stream``
---------------------

There are two structs defined in ``xml_io_struct.h``:

.. code-block::  cpp
  :caption: The original structs that should be overloaded
  :linenos:

  //! Name overloads
  template <typename T>
  struct xml_io_stream_name {
    static constexpr std::string_view name = "<unknown>"sv;
  };

  template <typename T>
  struct xml_io_stream {
    constexpr static std::string_view type_name = xml_io_stream_name_v<T>;

    // Core IO
    static void write(std::ostream&,
                      const T&,
                      bofstream*       = nullptr,
                      std::string_view = ""sv)                = delete;
    static void read(std::istream&, T&, bifstream* = nullptr) = delete;

    // Binary streaming IO (optional)
    static void get(std::span<const T>, bofstream*) = delete;
    static void put(std::span<T>, bifstream*)       = delete;

    // Text streaming (optional)
    static void parse(std::span<T>, std::istream&) = delete;
  };

The first struct allows overloading the ``type_name``.
This is a way to give a unique name to templated types.
For any generic type, this defaults to an unknown name,
which should never find itself into an actual file.

.. note::
  ``xml_io_stream_name_v<T>`` is the same as ``xml_io_stream_name<T>::name``.

The second struct contains the information required to access
the XML IO.  For the generic type, this has an unknown name
with all of its static methods deleted.  The intent is that
we need to overload these methods for specific types.
Please see a later subsection on helper structs that exists
that will allow you to skip writing overloads for some subtypes.

.. tip::
  As a reminder, to explicitly overload this structure
  you write ``template <> struct<MyType> {/* code */};``.  To template overload it
  you write ``template <typename ...T> struct<MyTemplate<T...>> {/* code */};``.

  The former always take precedence to the latter, so you can overload a generic template
  and then still specialize it.  Be very careful about the latter overload.  If the
  template matches any other overloaded template, there will be a name clash and the
  compiler will produce very complicated errors.

The ``type_name``, the ``write`` and the ``read`` method
are required to be overloaded for every workspace group in ARTS.
There is a convenience concept defined, ``arts_xml_ioable``,
which is used throughout the XML IO interface methods to ensure
that the name and methods are declared.  This should hopefully
ensure that accidentally missing to overload the methods will
cause a visible error directly in the code you write rather than
deep inside some generic template function.  In the past, this
type of error was often a link-time error with no visible
side-effect to the developer.

.. note::
  You will still get a linker error if your declaration does not come with a compiled definition.

The other three methods are there to aid template chaining.
You do not need to define them.  They may be used in manual
code as well, of course, and if your type fulfills the criterias
to define them, it is useful to do so.

First, ``put`` and ``get``.  If an overload defines these,
it fulfills the ``xml_io_binary`` concept.
These are for data types where the
memory layout is known when the type is defined.  ``Numeric``
is one such example.  ``Numeric`` will always be bit-by-bit
copyable and replacable.  You can store its data in binary
form without any issues.  In fact, doing so often results in
much faster file IO.  Since ``Numeric`` is binary-representable,
any data that can be represented as a ``std::span<Numeric>`` is
also binary-representable.  Such types involve ``std::array<Numeric, N>``
and ``std::vector<Numeric>`` and ``Vector`` and ``Matrix`` and so on.
Of these four types,  ``std::array<Numeric, N>`` is in turn also
binary-representable, so if we define ``put`` and ``get``
for it (spoiler: we do),
a ``std::vector<std::array<Numeric, N>>`` is also known
to be binary-representable.
``Vector``, however, is not binary-representable.
Its size is only known at runtime.

The last overloadable method is ``parse``.  When it is overloaded,
the concept ``xml_io_parseable`` true.  This method is supposed to
take a stream where we know there will be a span of data to read.
Again using ``Numeric`` as an example, we know it is going to be a single
connected set of characters that must be parsed through
the lense of ``double_imanip()``
or similar lenses to be understood.  Thus we know that if we have a span of
``Numeric`` to be parsed, they will be laid out next to eachother with some
whitespace separation, e.g. ``"1 2 3 4 5"``.  So it can be read.
Again, we can chain parse ``std::array``, but not ``std::vector``
since the latter has a runtime size.

It is very useful to overload the above whenever possible.
They unlock contiguous IO optimizations that are incredibly
useful to speed up reading routines.  On the other hand,
if your type cannot be represented as a span, it might not
help as much to force the overload.

Templated overloads of ``xml_io_stream``
----------------------------------------

Several types have template overloads that should handle the most common
combinations of advanced types in ARTS.  You are free to do an explicit
overload to any of these types, but are not required to.  The base classes
covers most use-cases.

These template overloads make use of the naming structure, and some types
have been explicitly renamed.  Only the base-name will be given below.
Overloaded names are too many.

What follows below is a summary of the existing template overloads and any
special considerations that must be had when wanting to use them.

``std::vector<T>`` or ``Array<T>``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are the same type and overloaded in the ``xml_io_stream_array.h`` file.
They will by default be named ``"Array"``.
The only limitation is that the type ``T`` must be ``arts_xml_ioable``.
If the ``T`` is parseable or binary-compatible, these methods are used.
In turn, the ``Array<T>`` will not be parseable or binary-compatible.

Example types:

- :class:`~pyarts3.arts.ArrayOfIndex`
- :class:`~pyarts3.arts.ArrayOfString`

``std::array<T, N>``
^^^^^^^^^^^^^^^^^^^^

These are also overloaded in the ``xml_io_stream_array.h`` file.
They will be default have the exact same name as  the ``Array<T>`` (even if
the ``Array<T>`` overload its name).
However, this name is not used in the XML file by defaunlt, instead opting to use
no name at all.
If the ``T`` is parseable or binary-compatible, these methods are used.
In turn, the ``std::array<T, N>`` will be parseable and/or binary-compatible
to the same degree as ``T``.

.. warning::
  You must overload the ``xml_io_stream<std::array<T, N>>`` if ``std::array<T, N>`` is
  directly facing the user and holds either parseable or binary compatible ``T``.
  This is because it does not respect XML tagging - it assumes that it is OK to just put its data where it is at.

All ``arts_options``
^^^^^^^^^^^^^^^^^^^^

These are overloaded to their names in ``xml_io_stream_enum_option.h``.
You can fullfill the criterias for this overload manually - or by a
very unfortunate mistake - by fullfilling the ``xml_enum_option`` concept.
These option values are put as strings in the XML-tag and parsed as such.

``std::function<T(Ts...)>``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are overloaded to the name ``"Function"`` in
``xml_io_stream_functional.h``.
Any attempt to call these overloads will result in
a runtime error.  The overload
exist so we can compose other XML IO in case there is a method in the future to
allow custom functional objects to be re-created.
Consider this compositional for all intent and purposes.

``std::unordered_map<Key, T>``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are overloaded to the name ``"Map"`` in ``xml_io_stream_map.h``.
Both ``Key`` and ``T`` must be ``arts_xml_ioable``.
They will be written one after the other for all types.
No care for binary-compatibility or parseability are taken.

``std::optional<T>``
^^^^^^^^^^^^^^^^^^^^

These are overloaded to the name ``"Optional"`` in
``xml_io_stream_optional.h``.
The type ``T`` must not only  be ``arts_xml_ioable``,
but also default constructible.
No care for binary-compatibility or parseability are taken.

``std::shared_ptr<T>``
^^^^^^^^^^^^^^^^^^^^^^

These are overloaded to the name ``"Shared"`` in
``xml_io_stream_shared_ptr.h``.
The non-constant type of ``T`` must not only  be
``arts_xml_ioable``, but also default constructible.
Constant types are first constructed as non-constant before
being moved into the constant shared pointer.
No care for binary-compatibility or parseability are taken.

``std::tuple<T...>``
^^^^^^^^^^^^^^^^^^^^

These are overloaded to the name ``"Tuple"`` in ``xml_io_stream_tuple.h``.
However, this name is not displayed to the user by default,
instead opting to use no name at all.
All types ``T...`` must  be ``arts_xml_ioable``.
Binary-compatibility and parseability is taken into account and forwarded
to nested types, with specialization if all ``T...`` are the same type,

.. warning::
  You must overload the ``xml_io_stream<std::tuple<T...>>`` if ``std::tuple<T...>`` is
  directly facing the user and if the tuple holds all parseable and binary compatible types.
  This is because it does not respect XML tagging - it assumes that it is OK to just put its data where it is at.

``std::pair<A, B>``
^^^^^^^^^^^^^^^^^^^

These are overloaded to do the same as
``std::tuple<A, B>``
in ``xml_io_stream_tuple.h``.
Both ``A`` and ``B`` must  be ``arts_xml_ioable``.
Binary-compatibility and parseability is taken into account and forwarded
to nested types, with specialization if both ``A`` and ``B`` are the same type,
if, and only if, both ``A`` and ``B`` share this property.

``std::variant<T...>``
^^^^^^^^^^^^^^^^^^^^^^

These are overloaded to the name of ``"Variant"`` in ``xml_io_stream_tuple.h``.
The types ``T...`` must not only be ``arts_xml_ioable``,
but also default constructible.
All the types ``T...`` must have unique ``xml_io_stream<T>::type_name...``.
These are used to identify which type is constructed
during the reading routine.
No care for binary-compatibility or parseability are taken.

Example types:

- :class:`~pyarts3.arts.AtmKeyVal`
- :class:`~pyarts3.arts.SurfKeyVal`

``matpack::data_t<T, N>``
^^^^^^^^^^^^^^^^^^^^^^^^^

These are overloaded to default name ``"Matpack"`` in
``xml_io_stream_matpack_mdspan.h``.
The type ``T`` must be ``arts_xml_ioable``.
Binary-compatibility and parseability is taken into account but not forwarded
to nested types.

Example types:

- :class:`~pyarts3.arts.Vector`
- :class:`~pyarts3.arts.Matrix`

``matpack::cdata_t<T, N...>``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are overloaded to the same name as the same-ranked and typed
``matpack::data_t<T, N>`` in ``xml_io_stream_matpack_mdspan.h``.
The type ``T`` must be ``arts_xml_ioable``.
Binary-compatibility and parseability is taken into account and forwarded
to nested types.

Example types:

- :class:`~pyarts3.arts.Vector2`
- :class:`~pyarts3.arts.Vector3`

``matpack::grid_t<Compare>``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Are completely forwarded as ``Vector`` in
``xml_io_stream_matpack_mdspan_helpers.h``.

Example types:

- :class:`~pyarts3.arts.AscendingGrid`
- :class:`~pyarts3.arts.DescendingGrid`

``matpack::gridded_data_t<T, Grids...>``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Are overloaded in ``xml_io_stream_matpack_mdspan_helpers.h``
to default name ``"GriddedField``.
The ``Grids...`` and the ``T`` must be ``arts_xml_ioable``.
No care for binary-compatibility or parseability are taken.

Example types:

- :class:`~pyarts3.arts.GriddedField1`
- :class:`~pyarts3.arts.GriddedField2`

Aggregates
^^^^^^^^^^

Are overloadable in ``xml_io_stream_aggregate.h``.  These will map directly
to the ``std::tuple<T...>`` overload when activated, except that it gives
type a named XML tag.

The overload is quite complicated, requiring somewhat hacky code since there
is currently no standard way to easily turn an aggregate into a tuple.
You activate it by modyfying the following piece of code for your type

.. code-block:: cpp
  :caption: The aggregate overload
  :linenos:

  // Example aggregate (do not copy)
  struct MyStruct {
    Index a, b;
  };

  // Overladable struct that activates the aggregate overlaods
  template<>
  struct xml_io_stream_aggregate<MyStruct> {
    static constexpr bool value = true;
  };

.. note::
  You will get the default naming scheme.  Please also overload the name.

.. warning::
  If any subtype of the aggregate is not possible to store, you will end
  up with very complicated errors.  Ensure all subtypes are XML IO compatible.

XML IO
======

This is the main way that we perform file IO in ARTS.

Each workspace group must support XML file IO.
If a group is added that does not support XML file IO, ARTS cannot compile.
Support is done by overloading a templated struct, ``xml_io_stream``,
that defines the core concepts involved in the workflow.

How the group handles XML file IO inside this is up to the group itself.
Generally, all lower level overloads should already be defined, so it
is common to just chain several calls to other XML IO methods.

The ``xml_io_stream``
---------------------

This is the core struct as defined in ``xml_io_struct.h`` as follows:

.. code-block::  cpp
  :caption: The original overload
  :linenos:

  template <typename T>
    struct xml_io_stream {
    constexpr static std::string_view type_name = xml_io_stream_name_v<T>;

    // Core IO
    static void write(std::ostream&,
                      const T&,
                      bofstream*       = nullptr,
                      std::string_view = ""sv)                = delete;
    static void read(std::istream&, T&, bifstream* = nullptr) = delete;

    // Binary/streaming IO (optional)
    static void get(const T* const, bofstream*, Size = 1) = delete;
    static void put(T*, bifstream*, Size = 1)             = delete;
    };

The struct explicitly deletes the overloads for the generic ``T`` type.
The file continues to define  the concepts ``arts_xml_ioable``, which
will be true for any overload ox ``xml_io_stream`` that crates
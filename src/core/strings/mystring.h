#pragma once

#include <array.h>
#include <configtypes.h>

#include <span>
#include <string>
#include <string_view>

#include "string_extract.h"

/** The String type for ARTS. Implementation. */
using String = std::string;

/** An array of Strings. */
using ArrayOfString = Array<String>;

/** An array of Strings. */
using ArrayOfArrayOfString = Array<Array<String>>;

/** Name string_view as we named string */
using StringView = std::string_view;

void tolower(String& x);

void toupper(String& x);

void join(String& res, const std::span<const String>& list, const String& with);
String join(const std::span<const String>& list, const String& with);

void split(ArrayOfString& aos, const String& x, const String& delim);
ArrayOfString split(const String& x, const String& delim);

void trim(String& x);

void replace(String& x, const String& from, const String& to);

//! Helper function when commas and spaces are needed after some first element
String comma(bool& first, const String& spaces = "");

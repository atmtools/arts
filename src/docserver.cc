/* Copyright (C) 2012 Oliver Lemke <olemke@core-dump.info>
 
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*!
  \file   docserver.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2010-09-21
 
  \brief  Implementation of the arts documentation server.
  */

#include "docserver.h"

#ifdef ENABLE_DOCSERVER

#include <stdint.h>
#include <algorithm>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "agenda_record.h"
#include "auto_version.h"
#include "global_data.h"
#include "libmicrohttpd/microhttpd.h"
#include "libmicrohttpd/platform.h"
#include "messages.h"
#include "methods.h"
#include "workspace_ng.h"

#define DOCSERVER_NAME "ARTS built-in documentation server"

#define DS_ERROR_404 "Page not found"

static int ahc_echo(void* cls,
                    struct MHD_Connection* connection,
                    const char* url,
                    const char* method,
                    const char* version,
                    const char* upload_data,
                    size_t* upload_data_size,
                    void** ptr);

//! Limit length of output
/*! Automatically inserts linebreaks at certain length.
 
  \author Oliver Lemke
  \date   2008-09-03
  */
void Docserver::limit_line_length(ostringstream& curline,
                                  ostringstream& token,
                                  const String& indent,
                                  size_t linelen) {
  if (indent.length() + curline.str().length() + token.str().length() >
      linelen) {
    get_os() << curline.str() << endl << indent;
    curline.str("");
  }
  curline << token.str();
  token.str("");
}

//! Intializes the docserver for a new page.
/** 
  Get the docserver ready for a new page.
 
  \param[in]      url          URL of page.
  \param[out]     String with content type of page.
 
  \author Oliver Lemke
  */
string Docserver::new_page(const string& url) {
  string surl = url;

  if (surl.find(get_baseurl()) == 0) surl.erase(0, get_baseurl().size());

  split_tokens(surl);

  while (tokens.size() && tokens[tokens.size() - 1] == "") tokens.pop_back();

  while (tokens.size() && tokens[0] == "") tokens.erase(tokens.begin());

  string content_type = "text/html; charset=utf-8";
  if (tokens.size() && tokens[tokens.size() - 1] == "styles.css") {
    insert_stylesheet();
    content_type = "text/css; charset=utf-8";
  } else if (tokens.size() && tokens[tokens.size() - 1] == "doccheck") {
    insert_broken_doc_links();
  } else {
    switch (tokens.size()) {
      case 0:
      case 1:
        insert_index();
        break;
      case 2:
        insert_doc();
        break;
      default:
        insert_error(DS_ERROR_404);
    }
  }

  return content_type;
}

//! Split url into tokens.
/** 
  Splits the string based on the given delimiter. The splitted string are
  appended to the elems vector.
 
  \param[in]      s          String to split.
 
  \author Oliver Lemke
  */
void Docserver::split_tokens(const string& s) {
  tokens.clear();

  stringstream ss(s);
  string item;
  while (getline(ss, item, '/')) tokens.push_back(item);
}

//! Convert character to HTML entity.
/** 
  Returns a string with either the original character or the HTML entity if
  it's a special char.
 
  \param[in]  ch     Character.
 
  \returns String with original character or HTML entity.
 
  \author Oliver Lemke
  */
string Docserver::html_escape_char(const char ch) {
  string ret;

  switch (ch) {
    case '<':
      ret.append("&lt;");
      break;
    case '>':
      ret.append("&gt;");
      break;
    default:
      ret.append(1, ch);
  }

  return ret;
}

//! Convert special characters in string to HTML entity.
/** 
  Returns a string with special characters replaced by their HTML entities.
 
  \param[in]  s     String.
 
  \returns String with original character or HTML entity.
 
  \author Oliver Lemke
  */
string Docserver::html_escape_string(const string& s) {
  string ret;

  for (string::const_iterator it = s.begin(); it != s.end(); it++)
    ret.append(html_escape_char(*it));

  return ret;
}

//! Open the content div.
/** 
  Output the starting content div tag to stream.
 
  \author Oliver Lemke
  */
void Docserver::begin_content() {
  get_os() << "<div class=\"content\">" << endl;
}

//! Close the content div.
/** 
  Output the ending content div tag to stream.
 
  \author Oliver Lemke
  */
void Docserver::end_content() { get_os() << "</div>" << endl; }

//! Output the docserver HTML page header to stream.
/** 
  Output the HTML header and start of the body to the stream.
 
  \param[in]      title  Page title.
 
  \author Oliver Lemke
  */
void Docserver::begin_page(string title) {
  if (title.length()) title += " - ";

  title += DOCSERVER_NAME;

  get_os()
      << "<!DOCTYPE html>" << endl
      << "<html lang=\"en\">" << endl
      << "<head>" << endl
      << "<title>" << title << "</title>" << endl
      << "<meta charset=\"utf-8\">" << endl
      << "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\" />"
      << endl
      << "<link rel=\"stylesheet\" href=\"" << mbaseurl << "/styles.css\">"
      << endl
      << "</head>" << endl
      << "<body>" << endl;
}

//! Output the docserver HTML page footer to stream.
/** 
  Output the HTML footer and end of the body to the stream.
 
  \author Oliver Lemke
  */
void Docserver::end_page() {
  get_os() << "<div class=\"footer\">Page generated by " << ARTS_FULL_VERSION
           << " - <a href=\"" << mbaseurl
           << "/doccheck\">Check docs</a></div>\n"
           << endl;

  get_os() << "</body>" << endl << "</html>";
}

//! Create HTML link to an agenda.
/** 
  Returns a string with the HTML link to the agenda.
 
  \param[in]  aname     Agenda name.
 
  \returns String containing the HTML link.
 
  \author Oliver Lemke
  */
String Docserver::insert_agenda_link(const String& aname) {
  ostringstream ret;
  ret << "<a href=\"" << mbaseurl << "/agendas/" << aname << "\">" << aname
      << "</a>";
  return ret.str();
}

//! Create HTML link to a workspace group.
/** 
  Returns a string with the HTML link to the workspace group.
 
  \param[in]  gname     Group name.
 
  \returns String containing the HTML link.
 
  \author Oliver Lemke
  */
String Docserver::insert_group_link(const String& gname) {
  ostringstream ret;
  ret << "<a href=\"" << mbaseurl << "/groups/" << gname << "\">" << gname
      << "</a>";
  return ret.str();
}

//! Create HTML link to a workspace method.
/** 
  Returns a string with the HTML link to the workspace method.
 
  \param[in]  mname     Method name.
 
  \returns String containing the HTML link.
 
  \author Oliver Lemke
  */
String Docserver::insert_wsm_link(const String& mname) {
  ostringstream ret;
  ret << "<a href=\"" << mbaseurl << "/methods/" << mname << "\">" << mname
      << "</a>";
  return ret.str();
}

//! Create HTML link to a workspace variable.
/** 
  Returns a string with the HTML link to the workspace variable.
 
  \param[in]  vname     Variable name.
 
  \returns String containing the HTML link.
 
  \author Oliver Lemke
  */
String Docserver::insert_wsv_link(const String& vname) {
  ostringstream ret;

  // Find wsv id:
  map<String, Index>::const_iterator it = Workspace::WsvMap.find(vname);
  if (it != Workspace::WsvMap.end()) {
    if (is_agenda_group_id(Workspace::wsv_data[it->second].Group()))
      ret << "<a href=\"" << mbaseurl << "/agendas/" << vname << "\">" << vname
          << "</a>";
    else
      ret << "<a href=\"" << mbaseurl << "/variables/" << vname << "\">"
          << vname << "</a>";
  }

  return ret.str();
}

//! Output list of agendas to stream.
/** 
  Output an HTML list of all agenda to the stream.
 
  \author Oliver Lemke
  */
void Docserver::list_agendas() {
  Index i;

  get_os() << "<h2>Agendas</h2>" << endl;

  get_os() << "<div class=\"firstcol\">" << endl << "<ul>" << endl;

  Index hitcount = 0;
  for (i = 0; i < Workspace::wsv_data.nelem(); ++i)
    if (is_agenda_group_id(Workspace::wsv_data[i].Group())) hitcount++;

  Index hitcount2 = 0;
  for (i = 0; i < Workspace::wsv_data.nelem(); ++i) {
    if (is_agenda_group_id(Workspace::wsv_data[i].Group())) {
      get_os() << "<li>" << insert_agenda_link(Workspace::wsv_data[i].Name())
               << "</li>" << endl;
      hitcount2++;

      if (hitcount2 == hitcount / 2)
        get_os() << "</ul>" << endl
                 << "</div>" << endl
                 << "<div class=\"secondcol\">" << endl
                 << "<ul>" << endl;
    }
  }

  get_os() << "</ul>" << endl << "</div>" << endl;
}

//! Output list of workspace groups to stream.
/** 
  Output an HTML list of all workspace groups to the stream.
 
  \author Oliver Lemke
  */
void Docserver::list_groups() {
  using global_data::wsv_group_names;
  Index i;

  get_os() << "<h2>Workspace Groups</h2>" << endl;

  get_os() << "<div class=\"firstcol\">" << endl << "<ul>" << endl;
  for (i = 0; i < wsv_group_names.nelem(); ++i) {
    get_os() << "<li>" << insert_group_link(wsv_group_names[i]) << "</li>"
             << endl;

    if (i + 1 == wsv_group_names.nelem() / 2)
      get_os() << "</ul>" << endl
               << "</div>" << endl
               << "<div class=\"secondcol\">" << endl
               << "<ul>" << endl;
  }

  get_os() << "</ul>" << endl << "</div>" << endl;
}

//! Output list of workspace methods to stream.
/** 
  Output an HTML list of all workspace methods to the stream.
 
  \author Oliver Lemke
  */
void Docserver::list_methods() {
  using global_data::md_data_raw;
  Index i;

  get_os() << "<h2>Workspace Methods</h2>" << endl;

  get_os() << "<div class=\"firstcol\">" << endl << "<ul>" << endl;
  for (i = 0; i < md_data_raw.nelem(); ++i) {
    get_os() << "<li>" << insert_wsm_link(md_data_raw[i].Name()) << "</li>"
             << endl;

    if (i + 1 == md_data_raw.nelem() / 2)
      get_os() << "</ul>" << endl
               << "</div>" << endl
               << "<div class=\"secondcol\">" << endl
               << "<ul>" << endl;
  }

  get_os() << "</ul>" << endl << "</div>" << endl;
}

//! Output list of workspace variables to stream.
/** 
  Output an HTML list of all workspace variables to the stream.
 
  \author Oliver Lemke
  */
void Docserver::list_variables() {
  Index i;

  get_os() << "<h2>Workspace Variables</h2>" << endl;

  get_os() << "<div class=\"firstcol\">" << endl << "<ul>" << endl;
  Index hitcount = 0;
  for (i = 0; i < Workspace::wsv_data.nelem(); ++i) {
    if (!is_agenda_group_id(Workspace::wsv_data[i].Group())) hitcount++;
  }

  Index hitcount2 = 0;
  for (i = 0; i < Workspace::wsv_data.nelem(); ++i) {
    if (!is_agenda_group_id(Workspace::wsv_data[i].Group())) {
      get_os() << "<li>" << insert_wsv_link(Workspace::wsv_data[i].Name())
               << "</li>" << endl;
      hitcount2++;

      if (hitcount2 == hitcount / 2)
        get_os() << "</ul>" << endl
                 << "</div>" << endl
                 << "<div class=\"secondcol\">" << endl
                 << "<ul>" << endl;
    }
  }

  get_os() << "</ul>" << endl << "</div>" << endl;
}

//! Substitute Workspace entities with links.
/** 
  Replaces WSVs, WSMs and groups enclosed in ** with html links.
 
  \param[in]      desc   String with the description.
  \param[in]      mname  Optional name of corresponding method.
 
  \returns        Description including HTML links.
 
  \author Oliver Lemke
  */
String Docserver::description_add_links(const String& desc,
                                        const String& mname) {
  string ret;
  string link;
  bool inside_link = false;
  string::const_iterator it = desc.begin();

  using global_data::AgendaMap;
  using global_data::MdRawMap;

  while (it != desc.end()) {
    if (!inside_link) {
      if (*it == '*')
        inside_link = true;
      else
        ret += html_escape_char(*it);
    } else {
      if (*it == '*') {
        inside_link = false;
        if (MdRawMap.find(link) != MdRawMap.end())
          ret += insert_wsm_link(link);
        else if (AgendaMap.find(link) != AgendaMap.end())
          ret += insert_agenda_link(link);
        else if (Workspace::WsvMap.find(link) != Workspace::WsvMap.end())
          ret += insert_wsv_link(link);
        else if (get_wsv_group_id(link) != -1)
          ret += insert_group_link(link);
        else if (mname != "") {
          using global_data::MdRawMap;
          bool found = false;

          // Find method id:
          map<String, Index>::const_iterator mit = MdRawMap.find(mname);
          if (mit != MdRawMap.end()) {
            using global_data::md_data_raw;
            const MdRecord& mdr = md_data_raw[mit->second];

            for (ArrayOfString::const_iterator sit = mdr.GIn().begin();
                 !found && sit != mdr.GIn().end();
                 sit++) {
              if ((*sit) == link) {
                ret += "*" + link + "*";
                found = true;
              }
            }

            for (ArrayOfString::const_iterator sit = mdr.GOut().begin();
                 !found && sit != mdr.GOut().end();
                 sit++) {
              if ((*sit) == link) {
                ret += "*" + link + "*";
                found = true;
              }
            }
          }

          if (!found)
            ret += "<span class=\"brokendoclink\">*" + link + "*</span>";
        } else
          ret += "<span class=\"brokendoclink\">*" + link + "*</span>";

        link = "";
      } else {
        if (!isalnum(*it) && *it != '_') {
          inside_link = false;
          ret += "*" + link + *it;
          link = "";
        } else
          link += html_escape_char(*it);
      }
    }

    it++;
  }

  if (inside_link) ret += "*" + link;

  return ret;
}

//! Output the workspace method documentation to stream.
/** 
  Output the documentation for the given workspace method to the stream in
  HTML formatting.
 
  \param[in]      mname  Method name.
 
  \author Oliver Lemke
  */
void Docserver::doc_method(const string& mname) {
  // Make global data visible:
  using global_data::md_data_raw;
  using global_data::MdRawMap;
  using global_data::wsv_group_names;

  // Let's first assume it is a method that the user wants to have
  // described.

  // Find method id:
  map<String, Index>::const_iterator it = MdRawMap.find(mname);
  if (it != MdRawMap.end()) {
    // If we are here, then the given name matches a method.
    const MdRecord& mdr = md_data_raw[it->second];
    String indent = "";

    get_os() << "<h3>Description</h3>" << endl;

    get_os() << "<pre>" << endl;
    get_os() << description_add_links(mdr.Description(), mname);
    get_os() << endl << "</pre>" << endl << endl;

    bool is_first_author = true;
    for (Index i = 0; i < mdr.Authors().nelem(); i++) {
      if (is_first_author) {
        get_os() << "<p><b>Authors: </b>";
        is_first_author = false;
      } else
        get_os() << ", ";

      get_os() << mdr.Authors()[i];
    }
    get_os() << "\n";

    // Print the method's synopsis
    while (indent.length() < mdr.Name().length() + 2) indent += ' ';

    get_os() << "<h3>Synopsis</h3>" << endl;

    ostringstream buf;
    ostringstream param;
    const size_t linelen = 2048;

    buf << "<p><table><tr><td>" << mdr.Name() << "(&nbsp;</td><td>";
    bool first = true;
    for (Index i = 0; i < mdr.Out().nelem(); ++i) {
      if (first)
        first = false;
      else
        buf << ", ";
      param << insert_wsv_link(Workspace::wsv_data[mdr.Out()[i]].Name());

      limit_line_length(buf, param, indent, linelen);
    }

    for (Index i = 0; i < mdr.GOutType().nelem(); ++i) {
      if (first)
        first = false;
      else
        buf << ", ";

      if (mdr.GOut()[i].length())
        param << mdr.GOut()[i];
      else
        param << "gout" << i;

      limit_line_length(buf, param, indent, linelen);
    }

    const ArrayOfIndex& inonly = mdr.InOnly();
    for (Index i = 0; i < inonly.nelem(); ++i) {
      if (first)
        first = false;
      else
        buf << ", ";

      param << insert_wsv_link(Workspace::wsv_data[inonly[i]].Name());

      limit_line_length(buf, param, indent, linelen);
    }

    for (Index i = 0; i < mdr.GInType().nelem(); ++i) {
      if (first)
        first = false;
      else
        buf << ", ";

      if (mdr.GIn()[i].length())
        param << mdr.GIn()[i];
      else
        param << "gin" << i;

      limit_line_length(buf, param, indent, linelen);
    }
    if (buf.str().length()) get_os() << buf.str();

    get_os() << " )</td></tr></table>" << endl;

    get_os() << "<h3>Variables</h3>" << endl;

    // Out:
    indent = "";
    String desc;
    get_os() << "<table>" << endl;
    for (Index i = 0; i < mdr.Out().nelem(); ++i) {
      buf.str("");
      buf << "<tr>";
      if (std::find(mdr.In().begin(), mdr.In().end(), mdr.Out()[i]) ==
          mdr.In().end())
        buf << "<td>OUT</td>";
      else
        buf << "<td>OUT+IN</td>";

      {
        const String& vname = Workspace::wsv_data[mdr.Out()[i]].Name();
        buf << "<td class=\"right\">" << insert_wsv_link(vname) << "</td><td>(";
        buf << insert_group_link(
            wsv_group_names[Workspace::wsv_data[mdr.Out()[i]].Group()]);
        buf << ")</td><td>";
      }

      get_short_wsv_description(
          desc, Workspace::wsv_data[mdr.Out()[i]].Description());

      if (buf.str().length() + desc.length() > linelen) {
        format_paragraph(desc, indent, linelen);
        buf << endl << indent << description_add_links(desc);
      } else {
        buf << description_add_links(desc);
      }

      get_os() << buf.str() << "</td></tr>" << endl;
    }

    size_t lastlen;
    bool fit;
    for (Index i = 0; i < mdr.GOut().nelem(); ++i) {
      buf.str("");
      buf << "<tr>";
      buf << "<td>GOUT</td><td class=\"right\">" << mdr.GOut()[i]
          << "</td><td>(";

      if (mdr.GOutType()[i] == get_wsv_group_id("Any") &&
          mdr.GOutSpecType()[i].nelem()) {
        bool firstarg = true;
        for (Index j = 0; j < mdr.GOutSpecType()[i].nelem(); j++) {
          if (!firstarg)
            buf << ", ";
          else
            firstarg = false;

          buf << insert_group_link(wsv_group_names[mdr.GOutSpecType()[i][j]]);
        }
      } else {
        buf << insert_group_link(wsv_group_names[mdr.GOutType()[i]]);
      }

      buf << ")</td><td>";
      desc = buf.str();
      lastlen = desc.length();
      fit = format_paragraph(desc, indent, linelen);
      buf.str("");
      get_os() << desc;

      desc = mdr.GOutDescription()[i];
      if (!fit) {
        format_paragraph(desc, indent, linelen);
        buf << endl << indent << description_add_links(desc);
      } else {
        if (lastlen + desc.length() > linelen) {
          format_paragraph(desc, indent, linelen, lastlen);
          buf << endl << description_add_links(desc);
        } else
          buf << description_add_links(desc);
      }

      get_os() << buf.str() << "</td></tr>" << endl;
    }

    for (Index i = 0; i < mdr.In().nelem(); ++i) {
      if (std::find(mdr.Out().begin(), mdr.Out().end(), mdr.In()[i]) !=
          mdr.Out().end())
        continue;

      buf.str("");
      buf << "<tr>";
      buf << "<td>IN</td>";

      const String& vname = Workspace::wsv_data[mdr.In()[i]].Name();
      buf << "<td class=\"right\">" << insert_wsv_link(vname);
      buf << "</td><td>(";
      buf << insert_group_link(
          wsv_group_names[Workspace::wsv_data[mdr.In()[i]].Group()]);
      buf << ")</td><td>";

      get_short_wsv_description(desc,
                                Workspace::wsv_data[mdr.In()[i]].Description());

      if (buf.str().length() + desc.length() > linelen) {
        format_paragraph(desc, indent, linelen, indent.length());
        buf << endl << indent << description_add_links(desc);
      } else
        buf << description_add_links(desc);

      get_os() << buf.str() << "</td></tr>" << endl;
    }

    for (Index i = 0; i < mdr.GIn().nelem(); ++i) {
      buf.str("");
      buf << "<tr>";
      buf << "<td>GIN</td><td class=\"right\">" << mdr.GIn()[i] << "</td><td>(";
      if (mdr.GInType()[i] == get_wsv_group_id("Any") &&
          mdr.GInSpecType()[i].nelem()) {
        bool firstarg = true;
        for (Index j = 0; j < mdr.GInSpecType()[i].nelem(); j++) {
          if (!firstarg)
            buf << ", ";
          else
            firstarg = false;
          buf << insert_group_link(wsv_group_names[mdr.GInSpecType()[i][j]]);
        }
      } else {
        buf << insert_group_link(wsv_group_names[mdr.GInType()[i]]);
      }

      if (mdr.GInDefault()[i] != NODEF) {
        buf << ", Default: ";
        if (mdr.GInType()[i] == get_wsv_group_id("String")) {
          buf << "\"" << mdr.GInDefault()[i] << "\"";
        } else {
          buf << mdr.GInDefault()[i];
        }
      }

      buf << ")</td><td>";
      desc = buf.str();
      lastlen = desc.length();
      fit = format_paragraph(desc, indent, linelen);
      buf.str("");
      get_os() << desc;

      desc = mdr.GInDescription()[i];
      if (!fit) {
        format_paragraph(desc, indent, linelen);
        buf << indent << description_add_links(desc);
      } else if (lastlen + desc.length() > linelen) {
        format_paragraph(desc, indent, linelen, indent.length());
        buf << indent << description_add_links(desc);
      } else {
        buf << description_add_links(desc);
      }

      get_os() << buf.str() << "</td></tr>" << endl;
    }
    get_os() << "</table>" << endl;
  } else {
    insert_error_message("There is no method by this name.");
  }
}

//! Output a list of methods that can use or generate the given variable.
/** 
  Output a list of all workspace methods and agendas that use, generate or
  require the given variable.
 
  \param[in]      vname  Variable name.
 
  \author Oliver Lemke
  */
void Docserver::doc_variable_methods(const string& vname) {
  // Check if the user gave the name of a specific variable.
  map<String, Index>::const_iterator mi = Workspace::WsvMap.find(vname);
  using global_data::md_data_raw;
  if (mi != Workspace::WsvMap.end()) {
    // If we are here, then the given name matches a variable.
    Index wsv_key = mi->second;
    Index hitcount = 0;

    // List specific methods:
    hitcount = 0;
    get_os() << "<h3>Specific methods that can generate " << vname << "</h3>"
             << endl
             << "<ul>" << endl;
    for (Index i = 0; i < md_data_raw.nelem(); ++i) {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];

      // This if statement checks whether Output, the list
      // of output variables contains the workspace
      // variable key.
      if (count(mdd.Out().begin(), mdd.Out().end(), wsv_key)) {
        get_os() << "<li>" << insert_wsm_link(mdd.Name()) << "\n";
        ++hitcount;
      }
    }
    if (0 == hitcount) get_os() << "<li>none\n";

    get_os() << endl << "</ul>" << endl;

    // List generic methods:
    get_os() << "<h3>Generic and supergeneric methods that can generate "
             << vname << "</h3>" << endl;
    get_os() << "<ul>" << endl;
    for (Index i = 0; i < md_data_raw.nelem(); ++i) {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];

      // This if statement checks whether GOutType, the list
      // of output variable types contains the group of the
      // requested variable.
      // The else clause picks up methods with supergeneric input.
      if (count(mdd.GOutType().begin(),
                mdd.GOutType().end(),
                Workspace::wsv_data[wsv_key].Group())) {
        get_os() << "<li>" << insert_wsm_link(mdd.Name()) << endl;
        ++hitcount;
      } else if (count(mdd.GOutType().begin(),
                       mdd.GOutType().end(),
                       get_wsv_group_id("Any"))) {
        for (Index j = 0; j < mdd.GOutType().nelem(); j++) {
          if (mdd.GOutType()[j] == get_wsv_group_id("Any")) {
            if (mdd.GOutSpecType()[j].nelem()) {
              if (count(mdd.GOutSpecType()[j].begin(),
                        mdd.GOutSpecType()[j].end(),
                        Workspace::wsv_data[wsv_key].Group())) {
                get_os() << "<li>" << insert_wsm_link(mdd.Name()) << endl;
                ++hitcount;
              }
            } else {
              get_os() << "<li>" << insert_wsm_link(mdd.Name()) << endl;
              ++hitcount;
            }
          }
        }
      }
    }
    if (0 == hitcount) get_os() << "<li>none" << endl;

    get_os() << endl << "</ul>" << endl;

    // List specific methods:
    hitcount = 0;
    get_os() << "<h3>Specific methods that require " << vname << "</h3>" << endl
             << "<ul>" << endl;
    for (Index i = 0; i < md_data_raw.nelem(); ++i) {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];

      // This if statement checks whether Output, the list
      // of output variables contains the workspace
      // variable key.
      if (count(mdd.In().begin(), mdd.In().end(), wsv_key)) {
        get_os() << "<li>" << insert_wsm_link(mdd.Name()) << "\n";
        ++hitcount;
      }
    }
    if (0 == hitcount) get_os() << "<li>none\n";

    get_os() << endl << "</ul>" << endl;

    // List generic methods:
    hitcount = 0;
    get_os() << "<h3>Generic and supergeneric methods that can use " << vname
             << "</h3>" << endl;
    get_os() << "<ul>" << endl;
    for (Index i = 0; i < md_data_raw.nelem(); ++i) {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];

      // This if statement checks whether GOutType, the list
      // of output variable types contains the group of the
      // requested variable.
      // The else clause picks up methods with supergeneric input.
      if (count(mdd.GInType().begin(),
                mdd.GInType().end(),
                Workspace::wsv_data[wsv_key].Group())) {
        get_os() << "<li>" << insert_wsm_link(mdd.Name()) << endl;
        ++hitcount;
      } else if (count(mdd.GInType().begin(),
                       mdd.GInType().end(),
                       get_wsv_group_id("Any"))) {
        for (Index j = 0; j < mdd.GInType().nelem(); j++) {
          if (mdd.GInType()[j] == get_wsv_group_id("Any")) {
            if (mdd.GInSpecType()[j].nelem()) {
              if (count(mdd.GInSpecType()[j].begin(),
                        mdd.GInSpecType()[j].end(),
                        Workspace::wsv_data[wsv_key].Group())) {
                get_os() << "<li>" << insert_wsm_link(mdd.Name()) << endl;
                ++hitcount;
              }
            } else {
              get_os() << "<li>" << insert_wsm_link(mdd.Name()) << endl;
              ++hitcount;
            }
          }
        }
      }
    }
    if (0 == hitcount) get_os() << "<li>none" << endl;

    get_os() << endl << "</ul>" << endl;

    // List agendas with this variable as output:
    using global_data::agenda_data;
    hitcount = 0;
    get_os() << "<h3>Agendas that can generate " << vname << "</h3>" << endl
             << "<ul>" << endl;
    for (Index i = 0; i < agenda_data.nelem(); ++i) {
      // Get handle on method record:
      const AgRecord& ar = agenda_data[i];

      // This if statement checks whether Output, the list
      // of output variables contains the workspace
      // variable key.
      if (count(ar.Out().begin(), ar.Out().end(), wsv_key)) {
        get_os() << "<li>" << insert_agenda_link(ar.Name()) << "\n";
        ++hitcount;
      }
    }
    if (0 == hitcount) get_os() << "<li>none\n";

    get_os() << endl << "</ul>" << endl;

    // List agendas with this variable as input:
    hitcount = 0;
    get_os() << "<h3>Agendas that require " << vname << "</h3>" << endl
             << "<ul>" << endl;
    for (Index i = 0; i < agenda_data.nelem(); ++i) {
      // Get handle on method record:
      const AgRecord& ar = agenda_data[i];

      // This if statement checks whether Output, the list
      // of output variables contains the workspace
      // variable key.
      if (count(ar.In().begin(), ar.In().end(), wsv_key)) {
        get_os() << "<li>" << insert_agenda_link(ar.Name()) << "\n";
        ++hitcount;
      }
    }

    if (0 == hitcount) get_os() << "<li>none\n";

    get_os() << endl << "</ul>" << endl;
  }
}

//! Output the workspace variable documentation to stream.
/** 
  Output the documentation for the given workspace variable to the stream in
  HTML formatting.
 
  \param[in]      vname  Variable name.
 
  \author Oliver Lemke
  */
void Docserver::doc_variable(const string& vname) {
  using global_data::wsv_group_names;

  // Find wsv id:
  map<String, Index>::const_iterator it = Workspace::WsvMap.find(vname);
  if (it != Workspace::WsvMap.end()) {
    // If we are here, then the given name matches a workspace
    // variable.
    get_os() << "<pre>" << endl;
    get_os() << description_add_links(
        Workspace::wsv_data[it->second].Description());
    get_os() << endl << "</pre>" << endl << endl;

    get_os() << "<p><b>Group: </b>"
             << insert_group_link(
                    wsv_group_names[Workspace::wsv_data[it->second].Group()])
             << endl;

    doc_variable_methods(vname);
  } else {
    insert_error_message("There is no variable by this name.");
  }
}

//! Output the agenda documentation to stream.
/** 
  Output the documentation for the given agenda to the stream in
  HTML formatting.
 
  \param[in]      aname  Agenda name.
 
  \author Oliver Lemke
  */
void Docserver::doc_agenda(const string& aname) {
  using global_data::wsv_group_names;

  // Find wsv id:
  map<String, Index>::const_iterator it = Workspace::WsvMap.find(aname);
  using global_data::agenda_data;
  using global_data::AgendaMap;
  map<String, Index>::const_iterator ait = AgendaMap.find(aname);

  if (it != Workspace::WsvMap.end() && ait != AgendaMap.end()) {
    // If we are here, then the given name matches a workspace
    // variable.
    get_os() << "<pre>" << endl;
    get_os() << description_add_links(agenda_data[ait->second].Description());
    get_os() << endl << "</pre>" << endl << endl;

    get_os() << "<p><b>Group: </b>"
             << insert_group_link(
                    wsv_group_names[Workspace::wsv_data[it->second].Group()])
             << endl;

    get_os() << "<h3>Variables</h3>" << endl;

    // Out:
    if (ait != AgendaMap.end()) {
      // If we are here, then the given name matches a method.
      const AgRecord& agr = agenda_data[ait->second];
      String indent = "";
      String desc;
      ostringstream buf;
      size_t linelen = 80;
      get_os() << "<table>" << endl;
      for (Index i = 0; i < agr.Out().nelem(); ++i) {
        buf.str("");
        buf << "<tr>";
        buf << "<td>OUT</td>";

        {
          const String& vname = Workspace::wsv_data[agr.Out()[i]].Name();
          buf << "<td class=\"right\">" << insert_wsv_link(vname)
              << "</td><td>(";
          buf << insert_group_link(
              wsv_group_names[Workspace::wsv_data[agr.Out()[i]].Group()]);
          buf << ")</td><td>";
        }

        get_short_wsv_description(
            desc, Workspace::wsv_data[agr.Out()[i]].Description());

        if (buf.str().length() + desc.length() > linelen) {
          format_paragraph(desc, indent, linelen);
          buf << endl << indent << description_add_links(desc);
        } else {
          buf << description_add_links(desc);
        }

        get_os() << buf.str() << "</td></tr>" << endl;
      }

      for (Index i = 0; i < agr.In().nelem(); ++i) {
        buf.str("");
        buf << "<tr>";
        buf << "<td>IN</td>";

        const String& vname = Workspace::wsv_data[agr.In()[i]].Name();
        buf << "<td class=\"right\">" << insert_wsv_link(vname);
        buf << "</td><td>(";
        buf << insert_group_link(
            wsv_group_names[Workspace::wsv_data[agr.In()[i]].Group()]);
        buf << ")</td><td>";

        get_short_wsv_description(
            desc, Workspace::wsv_data[agr.In()[i]].Description());

        if (buf.str().length() + desc.length() > linelen) {
          format_paragraph(desc, indent, linelen, indent.length());
          buf << endl << indent << description_add_links(desc);
        } else
          buf << description_add_links(desc);

        get_os() << buf.str() << "</td></tr>" << endl;
      }

      get_os() << "</table>" << endl;
    }

    doc_variable_methods(aname);
  } else {
    insert_error_message("There is no agenda by this name.");
  }
}

//! Output the workspace group documentation to stream.
/** 
  Output the documentation for the given workspace group to the stream in
  HTML formatting.
 
  \param[in]      gname  Group name.
 
  \author Oliver Lemke
  */
void Docserver::doc_group(const string& gname) {
  // Check if the user gave the name of a specific variable.
  Index gid = get_wsv_group_id(gname);
  using global_data::md_data_raw;
  if (gid != -1) {
    // If we are here, then the given name matches a group.
    Index hitcount = 0;

    if (gname != "Any") {
      // List specific methods:
      hitcount = 0;
      get_os() << "<h3>Specific methods that can generate " << gname << "</h3>"
               << endl;
      get_os() << "<ul>" << endl;
      for (Index i = 0; i < md_data_raw.nelem(); ++i) {
        // Get handle on method record:
        const MdRecord& mdd = md_data_raw[i];

        bool first = true;
        for (Index j = 0; j < mdd.Out().nelem(); j++) {
          // This if statement checks whether the type of this output variable
          // matches this group.
          if (Workspace::wsv_data[mdd.Out()[j]].Group() == gid) {
            if (first) {
              first = false;
              get_os() << "<li>" << insert_wsm_link(mdd.Name()) << " (";
            } else
              get_os() << ", ";
            get_os() << insert_wsv_link(
                Workspace::wsv_data[mdd.Out()[j]].Name());

            ++hitcount;
          }
        }
        if (!first) get_os() << ")" << endl;
      }
      if (0 == hitcount) get_os() << "<li>none" << endl;

      get_os() << endl << "</ul>" << endl;
    }

    // List generic methods:
    get_os() << "<h3>Generic and supergeneric methods that can generate "
             << gname << "</h3>" << endl;
    get_os() << "<ul>" << endl;
    for (Index i = 0; i < md_data_raw.nelem(); ++i) {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];

      // This if statement checks whether GOutType, the list
      // of output variable types contains the group of the
      // requested variable.
      // The else clause picks up methods with supergeneric input.
      if (count(mdd.GOutType().begin(), mdd.GOutType().end(), gid)) {
        get_os() << "<li>" << insert_wsm_link(mdd.Name()) << endl;
        ++hitcount;
      } else if (count(mdd.GOutType().begin(),
                       mdd.GOutType().end(),
                       get_wsv_group_id("Any"))) {
        for (Index j = 0; j < mdd.GOutType().nelem(); j++) {
          if (mdd.GOutType()[j] == get_wsv_group_id("Any")) {
            if (mdd.GOutSpecType()[j].nelem()) {
              if (count(mdd.GOutSpecType()[j].begin(),
                        mdd.GOutSpecType()[j].end(),
                        gid)) {
                get_os() << "<li>" << insert_wsm_link(mdd.Name()) << endl;
                ++hitcount;
              }
            } else {
              get_os() << "<li>" << insert_wsm_link(mdd.Name()) << endl;
              ++hitcount;
            }
          }
        }
      }
    }
    if (0 == hitcount) get_os() << "<li>none" << endl;

    get_os() << endl << "</ul>" << endl;

    if (gname != "Any") {
      hitcount = 0;
      get_os() << "<h3>Specific methods that require variables of group "
               << gname << "</h3>" << endl;
      get_os() << "<ul>" << endl;
      for (Index i = 0; i < md_data_raw.nelem(); ++i) {
        // Get handle on method record:
        const MdRecord& mdd = md_data_raw[i];

        bool first = true;
        for (Index j = 0; j < mdd.In().nelem(); j++) {
          // This if statement checks whether the type of this output variable
          // matches this group.
          if (Workspace::wsv_data[mdd.In()[j]].Group() == gid) {
            if (first) {
              first = false;
              get_os() << "<li>" << insert_wsm_link(mdd.Name()) << " (";
            } else
              get_os() << ", ";
            get_os() << insert_wsv_link(
                Workspace::wsv_data[mdd.In()[j]].Name());

            ++hitcount;
          }
        }
        if (!first) get_os() << ")" << endl;
      }
      if (0 == hitcount) get_os() << "<li>none" << endl;

      get_os() << endl << "</ul>" << endl;
    }

    hitcount = 0;
    get_os() << "<h3>Generic and supergeneric methods that can use " << gname
             << "</h3>" << endl;
    get_os() << "<ul>" << endl;
    for (Index i = 0; i < md_data_raw.nelem(); ++i) {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];

      // This if statement checks whether GOutType, the list
      // of output variable types contains the group of the
      // requested variable.
      // The else clause picks up methods with supergeneric input.
      if (count(mdd.GInType().begin(), mdd.GInType().end(), gid)) {
        get_os() << "<li>" << insert_wsm_link(mdd.Name()) << endl;
        ++hitcount;
      } else if (count(mdd.GInType().begin(),
                       mdd.GInType().end(),
                       get_wsv_group_id("Any"))) {
        for (Index j = 0; j < mdd.GInType().nelem(); j++) {
          if (mdd.GInType()[j] == get_wsv_group_id("Any")) {
            if (mdd.GInSpecType()[j].nelem()) {
              if (count(mdd.GInSpecType()[j].begin(),
                        mdd.GInSpecType()[j].end(),
                        gid)) {
                get_os() << "<li>" << insert_wsm_link(mdd.Name()) << endl;
                ++hitcount;
              }
            } else {
              get_os() << "<li>" << insert_wsm_link(mdd.Name()) << endl;
              ++hitcount;
            }
          }
        }
      }
    }
    if (0 == hitcount) get_os() << "<li>none" << endl;

    get_os() << endl << "</ul>" << endl;

    if (gname != "Any") {
      Index i;

      get_os() << "<h3>Workspace Variables of group " << gname << "</h3>"
               << endl
               << "<ul>" << endl;

      hitcount = 0;
      for (i = 0; i < Workspace::wsv_data.nelem(); ++i) {
        if (Workspace::wsv_data[i].Group() == get_wsv_group_id(gname)) {
          get_os() << "<li>" << insert_wsv_link(Workspace::wsv_data[i].Name())
                   << endl;
          hitcount++;
        }
      }
      if (0 == hitcount) get_os() << "<li>none" << endl;

      get_os() << "</ul>" << endl;
    }
  } else {
    insert_error_message("There is no group by this name.");
  }
}

//! Find token type.
/** 
  If the first element of Docserver::tokens is "all", look at the second token
  and try to determine whether it's a methods, variable, agenda or group. Then
  changes the first element of tokens accordingly.

  \author Oliver Lemke
  */
void Docserver::find_token_type() {
  if (tokens.size() < 1 || tokens[0] != "all") return;

  // Make global data visible:
  using global_data::AgendaMap;
  using global_data::MdRawMap;

  // Find method id:
  if (MdRawMap.find(tokens[1]) != MdRawMap.end())
    tokens[0] = "methods";
  else if (AgendaMap.find(tokens[1]) != AgendaMap.end())
    tokens[0] = "agendas";
  else if (Workspace::WsvMap.find(tokens[1]) != Workspace::WsvMap.end())
    tokens[0] = "variables";
  else if (get_wsv_group_id(tokens[1]) != -1)
    tokens[0] = "groups";
}

//! Output HTML code for a breadcrumb token.
/** 
  Output the HTML code for the token with the given id. Used for the
  navigation at the top of each docserver page.
 
  \param[in]      token_id   Id of the desired token.
 
  \author Oliver Lemke
  */
void Docserver::insert_breadcrumb_token(size_t token_id) {
  if (token_id != tokens.size()) {
    get_os() << "<a href=\"" << mbaseurl << "/";
    for (size_t t = 0; t < token_id; t++) {
      if (t) get_os() << "/";
      get_os() << tokens[t];
    }
    get_os() << "\">";
  }

  if (!token_id)
    get_os() << "Home";
  else if (tokens[token_id - 1] == "methods")
    get_os() << "Methods";
  else if (tokens[token_id - 1] == "variables")
    get_os() << "Variables";
  else if (tokens[token_id - 1] == "agendas")
    get_os() << "Agendas";
  else if (tokens[token_id - 1] == "groups")
    get_os() << "Groups";
  else if (tokens[token_id - 1] == "all")
    get_os() << "All";
  else if (tokens[token_id - 1] == "doccheck")
    get_os() << "Doc Check";
  else
    get_os() << tokens[token_id - 1];

  if (token_id != tokens.size()) get_os() << "</a>";
}

//! Output HTML code for breadcrumbs.
/** 
  Output the HTML code for breadcrumbs based on Docserver::tokens. Used for the
  navigation at the top of each docserver page.
 
  \author Oliver Lemke
  */
void Docserver::insert_breadcrumbs() {
  get_os() << "<div id=\"navbar\"><div class=\"breadcrumbs\">";
  for (size_t t = 0; t <= tokens.size(); t++) {
    if (t) get_os() << "&nbsp;>>&nbsp;";
    insert_breadcrumb_token(t);
  }
  get_os() << "</div>" << endl;

  get_os() << "<div class=\"goto\">Go to: "
           << "<a href=\"" << mbaseurl << "/groups/\">Groups</a>&nbsp;-&nbsp;"
           << "<a href=\"" << mbaseurl
           << "/variables/\">Variables</a>&nbsp;-&nbsp;"
           << "<a href=\"" << mbaseurl << "/methods/\">Methods</a>&nbsp;-&nbsp;"
           << "<a href=\"" << mbaseurl << "/agendas/\">Agendas</a>"
           << "</div></div>" << endl;
}

//! Output error.
/** 
  Outputs an error page.
 
  \param[in]      error   error string.
 
  \author Oliver Lemke
  */
void Docserver::insert_error(const string& error) {
  begin_page("");
  insert_breadcrumbs();
  begin_content();
  insert_error_message(error);
  end_content();
  end_page();
}

//! Output error message.
/** 
  Formats the given string as an error message.
 
  \param[in]      error  error string.
 
  \author Oliver Lemke
  */
void Docserver::insert_error_message(const string& error) {
  if (error.length())
    get_os() << "<p class=\"error\">" << error << "</p>" << endl;
}

//! Output H1 HTML tag.
/** 
  Formats the title string as an H1 header.
 
  \param[in]      title  Title string.
 
  \author Oliver Lemke
  */
void Docserver::insert_title(const string& title) {
  get_os() << "<h1>" DOCSERVER_NAME;
  if (title.length()) get_os() << " &mdash; " << title;
  get_os() << "</h1>" << endl;
}

//! Output an HTML index of the workspace members.
/** 
  Output depending on the tokens either a list of workspace methods, variables,
  agendas or groups. Or if tokens is empty, a list of all of them.
 
  \author Oliver Lemke
  */
void Docserver::insert_index() {
  if (tokens.size() == 0 || tokens[0] == "all") {
    begin_page("");
    insert_breadcrumbs();
    begin_content();
    insert_title("Index");

    list_groups();
    list_variables();
    list_methods();
    list_agendas();
    end_content();
    end_page();
  } else {
    if (tokens[0] == "methods") {
      begin_page("Method Index");
      insert_breadcrumbs();
      begin_content();
      insert_title("Method Index");
      get_os() << "<table class=\"list\">" << endl;
      list_methods();
      get_os() << "</table>" << endl;
      end_content();
      end_page();
    } else if (tokens[0] == "variables") {
      begin_page("Variable Index");
      insert_breadcrumbs();
      begin_content();
      insert_title("Variable Index");
      get_os() << "<table class=\"list\">" << endl;
      list_variables();
      get_os() << "</table>" << endl;
      end_content();
      end_page();
    } else if (tokens[0] == "groups") {
      begin_page("Group Index");
      insert_breadcrumbs();
      begin_content();
      insert_title("Group Index");
      get_os() << "<table class=\"list\">" << endl;
      list_groups();
      get_os() << "</table>" << endl;
      end_content();
      end_page();
    } else if (tokens[0] == "agendas") {
      begin_page("Agenda Index");
      insert_breadcrumbs();
      begin_content();
      insert_title("Agenda Index");
      get_os() << "<table class=\"list\">" << endl;
      list_agendas();
      get_os() << "</table>" << endl;
      end_content();
      end_page();
    } else {
      insert_error(DS_ERROR_404);
    }
  }
}

//! Output HTML documentation of a workspace member.
/**
  Output depending on the tokens the documentation of a workspace method,
  variable, agenda or group.

  \author Oliver Lemke
  */
void Docserver::insert_doc() {
  find_token_type();

  if (tokens[0] == "methods") {
    begin_page(tokens[1]);
    insert_breadcrumbs();
    begin_content();
    insert_title();
    get_os() << "<h2>"
             << "Workspace Method " + tokens[1] << "</h2>" << endl;
    doc_method(tokens[1]);
    end_content();
    end_page();
  } else if (tokens[0] == "variables") {
    begin_page(tokens[1]);
    insert_breadcrumbs();
    begin_content();
    insert_title();
    get_os() << "<h2>"
             << "Workspace Variable " + tokens[1] << "</h2>" << endl;
    doc_variable(tokens[1]);
    end_content();
    end_page();
  } else if (tokens[0] == "groups") {
    begin_page(tokens[1]);
    insert_breadcrumbs();
    begin_content();
    insert_title();
    get_os() << "<h2>"
             << "Workspace Group " + tokens[1] << "</h2>" << endl;
    doc_group(tokens[1]);
    end_content();
    end_page();
  } else if (tokens[0] == "agendas") {
    begin_page(tokens[1]);
    insert_breadcrumbs();
    begin_content();
    insert_title();
    get_os() << "<h2>"
             << "Agenda " + tokens[1] << "</h2>" << endl;
    doc_agenda(tokens[1]);
    end_content();
    end_page();
  } else {
    insert_error(DS_ERROR_404);
  }
}

//! Output docserver stylesheet.
/** 
  Outpus the docserver stylesheet to the stream.
 
  \author Oliver Lemke
  */
void Docserver::insert_stylesheet() {
  get_os()
      << "body { font-family: monospace; }" << endl
      << "a:link { color: #3465a4; text-decoration: none; }" << endl
      << "a:visited { color: #729fcf; text-decoration: none; }" << endl
      << "a:active { color: #ce5c00; text-decoration: none; background-color: #eeeeec}"
      << endl
      << "a:hover { color: #f57900; text-decoration: none; }" << endl

      << "@media (prefers-color-scheme: dark) {" << endl
      << "  body { background-color: #121212; color: #b0bec5; }" << endl
      << "  a:link { color: #90caf9; }" << endl
      << "  a:visited { color: #bbdefb; }" << endl
      << "  a:hover { color: #ffcc80; }" << endl
      << "  a:active { color: #ffcc80; background-color: #121212; }"
      << endl
      << "}" << endl

      << "table.list {" << endl
      << "width: 90%;" << endl
      << "margin-left: 5%;" << endl
      << "margin-right: 5%;" << endl
      << "}" << endl

      << "h1 {" << endl
      << "font-size: 1.5em;" << endl
      << "}" << endl

      << "h2 {" << endl
      << "font-size: 1.25em;" << endl
      << "}" << endl

      << "h3 {" << endl
      << "font-size: 1em;" << endl
      << "}" << endl

      << "li {" << endl
      << "font-size: 1em;" << endl
      << "}" << endl

      << "#navbar {" << endl
      << "position: fixed;" << endl
      << "top: 0px;" << endl
      << "left: 10px;" << endl
      << "right: 10px;" << endl
      << "background-color: #fff;" << endl
      << "border-bottom: solid 1px #ddd;" << endl
      << "border-left: solid 1px #ddd;" << endl
      << "border-right: solid 1px #ddd;" << endl
      << "padding: 2px;" << endl
      << "}" << endl

      << "@media (prefers-color-scheme: dark) {" << endl
      << "  #navbar {" << endl
      << "    background-color: #121212;" << endl
      << "    color: #b0bec5;" << endl
      << "    border-bottom-color: #b0bec5;" << endl
      << "    border-left-color: #b0bec5;" << endl
      << "    border-right-color: #b0bec5;" << endl
      << "  }" << endl
      << "}" << endl

      << ".firstcol {" << endl
      << "float: left;" << endl
      << "clear: left;" << endl
      << "width: 50%;" << endl
      << "white-space: nowrap;" << endl
      << "}" << endl

      << ".firstcol ul {" << endl
      << "float: left;" << endl
      << "clear: both;" << endl
      << "padding-top: 0;" << endl
      << "}" << endl

      << ".secondcol ul {" << endl
      << "float: left;" << endl
      << "clear: both;" << endl
      << "padding-top: 0;" << endl
      << "}" << endl

      << ".secondcol {" << endl
      << "float: left;" << endl
      << "clear: right;" << endl
      << "width: 50%;" << endl
      << "white-space: nowrap;" << endl
      << "}" << endl

      << ".firstcol ul li {" << endl
      << "margin-left: 0;" << endl
      << "}" << endl

      << ".brokendoclink {" << endl
      << "  color: #f44336;" << endl
      << "}" << endl

      << "@media (prefers-color-scheme: dark) {" << endl
      << "  .brokendoclink {" << endl
      << "    color: #ef9a9a;" << endl
      << "  }" << endl
      << "}" << endl

      << ".goto {" << endl
      << "font-size: small;" << endl
      << "float: right;" << endl
      << "}" << endl

      << ".breadcrumbs {" << endl
      << "font-size: small;" << endl
      << "float: left;" << endl
      << "}" << endl

      << "@media only screen and (max-device-width: 480px) {" << endl
      << "#navbar { position: static; border: none; }" << endl
      << ".goto { position: static; float: none; }" << endl
      << ".breadcrumbs { position: static; float: none; }" << endl
      << ".firstcol { float: left; clear: left; width: 100%; }" << endl
      << ".secondcol { float: left; clear: both; width: 100%; }" << endl
      << ".firstcol ul { margin-top: 0; margin-bottom: 0; }" << endl
      << ".secondcol ul { margin-top: 0; }" << endl
      << "ul { padding-left: 1em; }" << endl
      << "}" << endl

      << "table {" << endl
      << "border-width: 0px;" << endl
      << "}" << endl

      << "table td {" << endl
      << "vertical-align: top;" << endl
      << "}" << endl

      << "table td.right {" << endl
      << "text-align: right;" << endl
      << "}" << endl

      << ".content {" << endl
      << "padding-top: 1em;" << endl
      << "clear: both;" << endl
      << "width: 100%;" << endl
      << "}" << endl

      << ".error {" << endl
      << "color: #f44336;" << endl
      << "font-weight: bold;" << endl
      << "font-size: 1.2em;" << endl
      << "}" << endl

      << "@media (prefers-color-scheme: dark) {" << endl
      << "  .error {" << endl
      << "    color: #ef9a9a;" << endl
      << "  }" << endl
      << "}" << endl

      << "div.footer {" << endl
      << "float: left;" << endl
      << "text-align: right;" << endl
      << "color: #aaaaa8;" << endl
      << "font-size: small;" << endl
      << "clear: left;" << endl
      << "margin-top: 2em;" << endl
      << "width: 100%;" << endl
      << "}" << endl

      << endl;
}

//! Output list of broken links.
/** 
  Output a list of all ** links pointing to non-existent WSVs, WSMs, Agendas or
  groups.
 
  \author Oliver Lemke
  */
void Docserver::insert_broken_doc_links() {
  begin_page(tokens[1]);
  insert_breadcrumbs();
  begin_content();
  insert_title("Broken links");

  // Broken links in WSV descriptions
  bool first = true;
  for (Index i = 0; i < Workspace::wsv_data.nelem(); ++i) {
    vector<string> broken_links;
    broken_links =
        find_broken_description_links(Workspace::wsv_data[i].Description());

    if (broken_links.size()) {
      if (first) {
        first = false;
        get_os() << "<h2>Variable descriptions</h2>" << endl;
      }
      get_os() << "<p>" << insert_wsv_link(Workspace::wsv_data[i].Name())
               << ": ";
      for (vector<string>::iterator it = broken_links.begin();
           it != broken_links.end();
           it++) {
        if (it != broken_links.begin()) get_os() << ", ";
        get_os() << *it;
      }
      get_os() << "</p>" << endl;
    }
  }

  // Broken links in agenda descriptions
  using global_data::agenda_data;
  first = true;
  for (Array<AgRecord>::const_iterator ait = agenda_data.begin();
       ait != agenda_data.end();
       ait++) {
    vector<string> broken_links;
    broken_links = find_broken_description_links(ait->Description());

    if (broken_links.size()) {
      if (first) {
        first = false;
        get_os() << "<h2>Agenda descriptions</h2>" << endl;
      }
      get_os() << "<p>" << insert_agenda_link(ait->Name()) << ": ";
      for (vector<string>::iterator it = broken_links.begin();
           it != broken_links.end();
           it++) {
        if (it != broken_links.begin()) get_os() << ", ";
        get_os() << *it;
      }
      get_os() << "</p>" << endl;
    }
  }

  // Broken links in method descriptions
  using global_data::md_data_raw;

  first = true;
  for (Array<MdRecord>::const_iterator mit = md_data_raw.begin();
       mit != md_data_raw.end();
       mit++) {
    vector<string> broken_links;
    broken_links =
        find_broken_description_links(mit->Description(), mit->Name());

    if (broken_links.size()) {
      if (first) {
        first = false;
        get_os() << "<h2>Method descriptions</h2>" << endl;
      }
      get_os() << "<p>" << insert_wsm_link(mit->Name()) << ": ";
      for (vector<string>::iterator it = broken_links.begin();
           it != broken_links.end();
           it++) {
        if (it != broken_links.begin()) get_os() << ", ";
        get_os() << *it;
      }
      get_os() << "</p>" << endl;
    }
  }

  end_content();
  end_page();
}

//! Find broken links in documentation.
/** 
  Returns a list of broken ** links in Workspace documentation.
 
  \param[in]      desc   String with the description.
  \param[in]      mname  Optional name of corresponding method.
 
  \returns        Vector with broken links.
 
  \author Oliver Lemke
  */
vector<string> Docserver::find_broken_description_links(const String& desc,
                                                        const String& mname) {
  vector<string> broken_links;
  string ret;
  string link;
  bool inside_link = false;
  string::const_iterator it = desc.begin();

  using global_data::AgendaMap;
  using global_data::MdRawMap;

  while (it != desc.end()) {
    if (!inside_link) {
      if (*it == '*')
        inside_link = true;
      else
        ret += html_escape_char(*it);
    } else {
      if (*it == '*') {
        inside_link = false;
        bool found = false;
        if (MdRawMap.find(link) != MdRawMap.end() ||
            AgendaMap.find(link) != AgendaMap.end() ||
            Workspace::WsvMap.find(link) != Workspace::WsvMap.end() ||
            get_wsv_group_id(link) != -1)
          found = true;
        else if (mname != "") {
          using global_data::MdRawMap;

          // Find method id:
          map<String, Index>::const_iterator mit = MdRawMap.find(mname);
          if (mit != MdRawMap.end()) {
            using global_data::md_data_raw;
            const MdRecord& mdr = md_data_raw[mit->second];

            for (ArrayOfString::const_iterator sit = mdr.GIn().begin();
                 !found && sit != mdr.GIn().end();
                 sit++) {
              if ((*sit) == link) {
                ret += "*" + link + "*";
                found = true;
              }
            }

            for (ArrayOfString::const_iterator sit = mdr.GOut().begin();
                 !found && sit != mdr.GOut().end();
                 sit++) {
              if ((*sit) == link) {
                ret += "*" + link + "*";
                found = true;
              }
            }
          }
        }

        if (!found)
          broken_links.push_back("<span class=\"brokendoclink\">*" + link +
                                 "*</span>");

        link = "";
      } else {
        if (!isalnum(*it) && *it != '_') {
          inside_link = false;
          ret += "*" + link + *it;
          link = "";
        } else
          link += html_escape_char(*it);
      }
    }

    it++;
  }

  //  if (inside_link) ret += "*" + link;

  return broken_links;
}

//! Starts docserver.
/** 
  Starts the HTTP server daemon.
 
  \param[in]  daemon   Flag to run as daemon or in the foreground.
 
  \returns Status code.
 
  \author Oliver Lemke
  */
int Docserver::launch(bool daemon) {
  struct MHD_Daemon* d;

  d = MHD_start_daemon(MHD_USE_THREAD_PER_CONNECTION | MHD_USE_DEBUG,
                       (uint16_t)mport,
                       NULL,
                       NULL,
                       &ahc_echo,
                       (void*)this,
                       MHD_OPTION_END);

  if (d == NULL) {
    cerr << "Error: Cannot start server. Maybe port " << mport
         << " is already in use?\n";
    return 1;
  } else {
    if (daemon)
      cerr << "ARTS docserver listening at http://localhost:" << mport << "\n";
    else
      cerr << "\n"
           << "===========================================================\n"
           << "Now point your web browser to http://localhost:" << mport << "\n"
           << "===========================================================\n\n"
           << "Press enter to exit.\n";
  }

  if (daemon)
    pause();
  else {
    (void)getc(stdin);
    cout << "Stopping docserver.\n";
    MHD_stop_daemon(d);
    cout << "Goodbye.\n";
  }
  return 0;
}

//! Construct docserver.
/** 
  Constructs the HTTP server daemon object.
 
  \param[in]  port     Port to listen on.
  \param[in]  baseurl  URL to prepend to all links.
 
  \author Oliver Lemke
  */
Docserver::Docserver(const Index port, const string& baseurl) : mos(NULL) {
  mbaseurl = baseurl;
  if (port == -1)
    mport = 9000;
  else
    mport = port;
}

//! HTML request responder.
/** 
  Callback function to serve the HTML pages. Based on an example from the
  libmicrohttpd library.
 
  \param[in,out]   cls               Docserver object.
  \param[in,out]   connection        Connection info.
  \param[in]       url               Requested URL.
  \param[in]       method            Request method, must be GET.
  \param[in]       version           Unused parameter.
  \param[in]       upload_data       Unused parameter.
  \param[in]       upload_data_size  Unused parameter.
  \param[in,out]  ptr                Call state.

  \returns Status code.

  \author Oliver Lemke
  */
static int ahc_echo(void* cls,
                    struct MHD_Connection* connection,
                    const char* url,
                    const char* method,
                    const char* version _U_,
                    const char* upload_data _U_,
                    size_t* upload_data_size _U_,
                    void** ptr) {
  static int aptr;
  string surl(url);
  struct MHD_Response* response;
  int ret;

  if (!cls) {
    cerr << "Docserver error: Docserver object reference is NULL.\n";
    return MHD_NO; /* unexpected method */
  }

  // Make a local copy of the docserver object for thread-safety
  Docserver docserver = *((Docserver*)cls);

  if (0 != strcmp(method, "GET")) {
    cerr << "Docserver error: Unexpected method " << method << ".\n";
    return MHD_NO; /* unexpected method */
  }
  if (&aptr != *ptr) {
    /* do never respond on first call */
    *ptr = &aptr;
    return MHD_YES;
  }
  *ptr = NULL; /* reset when done */
  MHD_lookup_connection_value(connection, MHD_GET_ARGUMENT_KIND, "q");

  string content_type;
  ostringstream hout;
  docserver.set_ostream(hout);
  content_type = docserver.new_page(surl);
  docserver.clear_ostream();

  response = MHD_create_response_from_data(
      hout.str().length(), (void*)hout.str().c_str(), MHD_NO, MHD_YES);

  if (response == NULL) {
    cerr << "Docserver error: response = 0\n";
    return MHD_NO;
  }

  MHD_add_response_header(response, "Content-type", content_type.c_str());
  ret = MHD_queue_response(connection, MHD_HTTP_OK, response);
  MHD_destroy_response(response);

  return ret;
}

#endif /* ENABLE_DOCSERVER */

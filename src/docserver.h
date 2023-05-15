/*!
  \file   docserver.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2010-09-21

  \brief  Declarations for the arts documentation server.
*/

#ifndef docserver_h
#define docserver_h

#include "arts.h"
#include "exceptions.h"

#ifdef ENABLE_DOCSERVER
class Docserver {
 private:
  string mbaseurl;
  Index mport;
  ostream* mos;
  vector<string> tokens;

  void begin_page(string title);
  void end_page();
  void begin_content();
  void end_content();

  String insert_agenda_link(const String& aname);
  String insert_group_link(const String& gname);
  String insert_wsm_link(const String& mname);
  String insert_wsv_link(const String& vname);

  void insert_title(const string& title = "");
  void insert_breadcrumb_token(size_t token_id);
  void insert_breadcrumbs();
  void insert_error_message(const string& error = "");

  void insert_stylesheet();
  void insert_index();
  void insert_doc();
  void insert_error(const string& error);

  void list_agendas();
  void list_groups();
  void list_methods();
  void list_variables();

  String description_add_links(const String& desc, const String& mname = "");

  void doc_method(const string& mname);
  void doc_variable_methods(const string& vname);
  void doc_variable(const string& vname);
  void doc_agenda(const string& aname);
  void doc_group(const string& gname);

  void find_token_type();

  static string html_escape_char(const char ch);
  static string html_escape_string(const string& s);

  void split_tokens(const string& s);

  void limit_line_length(ostringstream& curline,
                         ostringstream& token,
                         const String& indent,
                         size_t linelen);

  ostream& get_os() {
    if (!mos) throw runtime_error("Output stream for docserver is NULL.");
    return *mos;
  }

 public:
  Docserver(const Index port, const string& baseurl = "");

  static std::tuple<size_t, std::vector<string> >
  list_broken_description_links();

  static vector<string> find_broken_description_links(const String& desc,
                                                      const String& mname = "");

  string new_page(const string& url);

  void set_ostream(ostream& os) { mos = &os; }
  void clear_ostream() { mos = NULL; }

  const string& get_baseurl() { return mbaseurl; }
  int launch(bool daemon);
};

#endif /* ENABLE_DOCSERVER */

void run_docserver(Index port = 9000,
                   const String& baseurl = "",
                   bool daemon = false);

#endif /* docserver_h */

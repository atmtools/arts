#include "arts.h"
#include "messages.h"

/** The output path. For example for the report file. */
string out_path;

/** The report file. */
ofstream report_file;

/** Verbosity levels. 
    @see Messages */
Messages messages;

//--------------------< The different output streams >--------------------

/** Level 0 output stream. @see OutStream */
Out0 out0;
/** Level 1 output stream. @see OutStream */
Out1 out1;
/** Level 2 output stream. @see OutStream */
Out2 out2;
/** Level 3 output stream. @see OutStream */
Out3 out3;

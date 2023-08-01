#include <algorithm>

#include "auto_wsg.h"
#include "workspace.h"

void append(Vector& v_out, const Vector& v_in, const Index& n) {
  append(v_out, v_in, static_cast<Numeric>(n));
}

void append(Vector& v_out, const Vector& v_in, const Numeric& n) {
  Vector v(v_in.size() + 1);
  std::copy(v_in.begin(), v_in.end(), v.begin());
  v.back() = n;
  v_out = v;
}

void add(Vector& y, const Numeric& x) { y += x; }

void new_y(Vector& y, const Numeric& x) { y = Vector{x}; }

template <>
void print<Index>(const Index&) {
  std::cout << "ind" << '\n';
}

template <>
void print<Any>(const Any&) {
  std::cout << "any" << '\n';
}

template <>
void print<Vector>(const Vector& v) {
  std::cout << v << '\n';
}

template <>
void print<Numeric>(const Numeric&) {
  std::cout << "Numeric" << '\n';
}

void create(Index&) {}

int main() {
  Workspace ws;
  //ws.set("y", Numeric{1.0});
  ws.set("x", Numeric{1.1});

  const Method printy("print", {}, {{"v", "y"}});
  const Method printx("print", {}, {{"v", "x"}});
  const Method printnum("print", {}, {{"v", "num"}});
  const Method new_y("new_y");
  const Method create("create", {"num"});
  const Method add("add", {"y", "x"});
  const Method append("append", {"y", "y"}, {{"n", "x"}});

  //add(ws);
  new_y(ws);
  std::cout << "y: " << ws.get<Vector>("y")
            << ";             0 x add(ws); 0 x append(ws)" << '\n';
  add(ws);
  std::cout << "y: " << ws.get<Vector>("y")
            << ";             1 x add(ws); 0 x append(ws)" << '\n';
  add(ws);
  std::cout << "y: " << ws.get<Vector>("y")
            << ";             2 x add(ws); 0 x append(ws)" << '\n';
  add(ws);
  std::cout << "y: " << ws.get<Vector>("y")
            << ";             3 x add(ws); 0 x append(ws)" << '\n';

  append(ws);
  std::cout << "y: " << ws.get<Vector>("y")
            << ";         3 x add(ws); 1 x append(ws)" << '\n';
  append(ws);
  std::cout << "y: " << ws.get<Vector>("y")
            << ";     3 x add(ws); 2 x append(ws)" << '\n';
  append(ws);
  std::cout << "y: " << ws.get<Vector>("y") << "; 3 x add(ws); 3 x append(ws)"
            << '\n';
  add(ws);
  std::cout << "y: " << ws.get<Vector>("y") << "; 4 x add(ws); 3 x append(ws)"
            << '\n';

  printy(ws);
  printx(ws);
  //printn(ws);
  create(ws);
  printnum(ws);
}
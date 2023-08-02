#include <algorithm>

#include "auto_wsa.h"
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
void print<Agenda>(const Agenda&) {
  std::cout << "Agenda" << '\n';
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
  const Method printy("print", {}, {{"v", "y"}});
  const Method printx("print", {}, {{"v", "x"}});
  const Method printnum("print", {}, {{"v", "num"}});
  const Method new_y("new_y");
  const Method create("create", {"num"});
  const Method add("add", {"y", "x"});
  const Method append("append", {"y", "y"}, {{"n", "x"}});
  const Method myy("y", Vector{1,2,3,4,5,6,7,8,9,10});


  Agenda y_out("y_out_agenda");
  y_out.add(myy);
  y_out.add(append);
  y_out.add(add);
  y_out.add(append);
  y_out.add(add);
  y_out.add(append);
  y_out.finalize();

/*
  {
    Workspace ws;
    //ws.set("y", Numeric{1.0});
    ws.set("x", Numeric{1.1});

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
  */
  
  {
    Workspace ws;
    y_out.execute(ws);
    printy(ws);
  }

  {
    Vector y{1,2,3,4};
    Workspace ws;
    y_out_agendaExecute(ws, y, y_out);
    std::cout << y << '\n';
  }
  
  {
    Workspace ws;
    y_out.execute(ws);
    printy(ws);
  }

  {
    Vector y{1,2,3,4};
    Workspace ws;
    y_out_agendaExecute(ws, y, y_out);
    std::cout << y << '\n';
  }

std::visit([](auto&& x){
  std::cout << *x << '\n';}, myy.get_setval().value().value);
}

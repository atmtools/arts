#include <workspace.h>

#include <cstdlib>
#include <iostream>
#include <optional>
#include <stdexcept>

using Input  = Generic<const Numeric, const Vector>;
using Output = Generic<Numeric, Vector>;
static_assert(std::variant_size_v<Input::Base> == 2);
static_assert(std::same_as<std::variant_alternative_t<0, Input::Base>, std::shared_ptr<const Numeric>>);
static_assert(std::same_as<std::variant_alternative_t<0, Output::Base>, std::shared_ptr<Numeric>>);
static_assert(std::variant_size_v<AnyInput::Base> > 100);
static_assert(std::is_convertible_v<const Wsv&, Input>);
static_assert(std::is_convertible_v<const Wsv&, Output>);
template <typename... Ts>
concept ValidGeneric = requires { typename Generic<Ts...>; };
template <typename T>
concept DeducibleGeneric = requires(T& value) { Generic(value); };
static_assert(not DeducibleGeneric<std::variant<Numeric, Vector>>);
static_assert(not DeducibleGeneric<std::variant<int>>);
static_assert(not ValidGeneric<int>);
static_assert(WorkspaceGroup<Numeric>);
static_assert(not WorkspaceGroup<const Numeric>);
static_assert(std::same_as<AnyOutput::Const, AnyInput>);
static_assert(ValidGeneric<Numeric, Vector>);
static_assert(ValidGeneric<const Numeric, const Vector>);
static_assert(std::constructible_from<Input, Numeric&>);
static_assert(std::constructible_from<Input, const Numeric&>);
static_assert(std::constructible_from<Input, Numeric>);
static_assert(std::constructible_from<Input, std::shared_ptr<Numeric>>);
static_assert(std::constructible_from<Input, std::shared_ptr<const Numeric>>);
static_assert(not std::constructible_from<Output, const Numeric&>);
static_assert(not std::constructible_from<Output, std::shared_ptr<const Numeric>>);
static_assert(std::same_as<decltype(*std::get<0>(std::declval<const Output&>())), Numeric&>);
static_assert(std::same_as<decltype(*std::get<0>(std::declval<const Input&>())), const Numeric&>);

static_assert(std::constructible_from<Output, const Generic<Numeric>&>);
static_assert(std::constructible_from<Input, const Output&>);
static_assert(not std::constructible_from<Output, const Input&>);
static_assert(not std::constructible_from<Generic<Numeric>, const Output&>);
static_assert(not std::constructible_from<Output, std::variant<Numeric, Vector>&>);
static_assert(not std::constructible_from<Input, std::variant<Numeric, Vector>&>);
static_assert(not std::constructible_from<Input, std::variant<Numeric, Vector>&&>);
static_assert(not std::constructible_from<Output, const std::variant<Numeric, Vector>&>);
static_assert(not std::constructible_from<Output, const std::variant<Numeric, Vector>&&>);

static_assert(std::same_as<Output::Const, Input>);
static_assert(std::same_as<Input::Const, Input>);

int main() try {
  std::optional<Input>   const_view;
  std::weak_ptr<Numeric> const_lifetime;
  {
    Output mutable_value(std::make_shared<Numeric>(12));
    const_lifetime = std::get<0>(mutable_value);
    if (std::get<0>(*const_view).get() != std::get<0>(mutable_value).get())
      throw std::runtime_error("Const conversion copied its pointee");
    *std::get<0>(mutable_value) = 13;
  }
  if (const_lifetime.expired() or *std::get<0>(*const_view) != 13)
    throw std::runtime_error("Const conversion lost ownership");
  const_view.reset();
  if (not const_lifetime.expired()) throw std::runtime_error("Const conversion leaked ownership");
  Generic<Numeric, Vector> mutable_source(std::make_shared<Numeric>(14));

  Output                owned(std::make_shared<Vector>(Vector{8, 9}));
  std::weak_ptr<Vector> retained = std::get<std::shared_ptr<Vector>>(owned);
  const Input           widened(std::move(owned));
  if (retained.expired() or std::get<std::shared_ptr<const Vector>>(widened).get() != retained.lock().get())
    throw std::runtime_error("Widening lost shared ownership");
  Numeric     borrowed  = 2;
  const Input read_only = borrowed;
  borrowed              = 3;
  if (*std::get<0>(read_only) != 3) throw std::runtime_error("Borrowed input copied its value");
  auto        vector       = std::make_shared<Vector>(Vector{1, 2});
  const Input const_vector = Wsv{std::shared_ptr<Vector>(vector)};
  if (std::get<1>(const_vector).get() != vector.get()) throw std::runtime_error("Const alternative lost identity");
  auto                 pointer = std::make_shared<Numeric>(4);
  std::optional<Input> input;
  {
    Wsv owner{std::shared_ptr<Numeric>(pointer)};
    input.emplace(owner);
    if (std::get<0>(*input).get() != pointer.get()) throw std::runtime_error("Input copied its value");
    Output output        = owner;
    *std::get<0>(output) = 7;
    if (*pointer != 7) throw std::runtime_error("Output lost workspace identity");
    AnyInput any = owner;
    if (std::get<std::shared_ptr<const Numeric>>(any).get() != pointer.get())
      throw std::runtime_error("Any input copied its value");
  }
  std::weak_ptr<Numeric> lifetime = pointer;
  pointer.reset();
  if (lifetime.expired() or *std::get<0>(*input) != 7) throw std::runtime_error("Input lost ownership");
  input.reset();
  if (not lifetime.expired()) throw std::runtime_error("Input leaked ownership");
  try {
    const Input invalid = Wsv{String{"wrong type"}};
    (void)invalid;
    throw std::logic_error("Unsupported alternative accepted");
  } catch (const std::runtime_error&) {}
  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return EXIT_FAILURE;
}

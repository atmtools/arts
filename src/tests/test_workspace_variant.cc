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
template <typename... Ts>
concept ValidGeneric = requires { typename Generic<Ts...>; };
static_assert(ValidGeneric<Numeric, Vector>);
static_assert(ValidGeneric<const Numeric, const Vector>);
static_assert(not ValidGeneric<Numeric, const Vector>);
static_assert(not ValidGeneric<const Numeric, Vector>);
static_assert(not ValidGeneric<Numeric, const Numeric>);
static_assert(std::constructible_from<Input, Numeric&>);
static_assert(std::constructible_from<Input, const Numeric&>);
static_assert(std::constructible_from<Input, Numeric>);
static_assert(std::constructible_from<Input, std::shared_ptr<Numeric>>);
static_assert(std::constructible_from<Input, std::shared_ptr<const Numeric>>);
static_assert(not std::constructible_from<Output, const Numeric&>);
static_assert(not std::constructible_from<Output, std::shared_ptr<const Numeric>>);
static_assert(std::same_as<decltype(*std::get<0>(std::declval<const Output&>())), Numeric&>);
static_assert(std::same_as<decltype(*std::get<0>(std::declval<const Input&>())), const Numeric&>);

using Extended = Generic<Numeric>::Extend<Generic<Numeric, Vector>, Generic<Vector, String>>;
static_assert(std::same_as<Extended, Generic<Numeric, Vector, String>>);
static_assert(std::constructible_from<Output, const Generic<Numeric>&>);
static_assert(std::constructible_from<Input, const Output&>);
static_assert(not std::constructible_from<Output, const Input&>);
static_assert(not std::constructible_from<Generic<Numeric>, const Output&>);
static_assert(not std::constructible_from<Output, const std::variant<Numeric, Vector>&>);
static_assert(not std::constructible_from<Output, const std::variant<Numeric, Vector>&&>);

static_assert(std::same_as<Output::Const, Input>);
static_assert(std::same_as<Input::Const, Input>);
static_assert(std::same_as<Input::Extend<Generic<const Vector, const String>>,
                           Generic<const Numeric, const Vector, const String>>);

int main() try {
  std::optional<Input>   const_view;
  std::weak_ptr<Numeric> const_lifetime;
  {
    Output mutable_value(std::make_shared<Numeric>(12));
    const_lifetime = std::get<0>(mutable_value);
    const_view     = mutable_value.as_const();
    if (std::get<0>(*const_view).get() != std::get<0>(mutable_value).get())
      throw std::runtime_error("Const conversion copied its pointee");
    *std::get<0>(mutable_value) = 13;
  }
  if (const_lifetime.expired() or *std::get<0>(*const_view) != 13)
    throw std::runtime_error("Const conversion lost ownership");
  const_view.reset();
  if (not const_lifetime.expired()) throw std::runtime_error("Const conversion leaked ownership");
  Generic<Numeric, Vector> mutable_source(std::make_shared<Numeric>(14));
  auto                     moved_const = std::move(mutable_source).as_const();
  if (std::get<0>(mutable_source) or *std::get<0>(moved_const) != 14)
    throw std::runtime_error("Const conversion did not transfer ownership");

  AtmKeyVal key = AtmKey::t;
  Generic   key_view(key);
  *std::get<std::shared_ptr<AtmKey>>(key_view) = AtmKey::p;
  if (std::get<AtmKey>(key) != AtmKey::p) throw std::runtime_error("Variant lvalue was copied");
  const AtmKeyVal water = SpeciesEnum::Water;
  Generic         water_view(water);
  static_assert(std::same_as<typename decltype(water_view)::Base,
                             Generic<const AtmKey,
                                     const SpeciesEnum,
                                     const SpeciesIsotope,
                                     const QuantumLevelIdentifier,
                                     const ScatteringSpeciesProperty>::Base>);
  if (std::get<std::shared_ptr<const SpeciesEnum>>(water_view).get() != &std::get<SpeciesEnum>(water))
    throw std::runtime_error("Const variant lvalue was copied");
  auto owned = Generic(std::variant<Numeric, Vector>{Vector{8, 9}});
  if ((*std::get<std::shared_ptr<Vector>>(owned))[1] != 9) throw std::runtime_error("Temporary variant lost its value");
  std::weak_ptr<Vector> retained = std::get<std::shared_ptr<Vector>>(owned);
  const Input           widened(std::move(owned));
  if (retained.expired() or std::get<std::shared_ptr<const Vector>>(widened).get() != retained.lock().get())
    throw std::runtime_error("Widening lost shared ownership");
  auto const_owned = Generic(std::move(water));
  if (std::get<std::shared_ptr<const SpeciesEnum>>(const_owned).get() == &std::get<SpeciesEnum>(water))
    throw std::runtime_error("Const rvalue variant was borrowed");
  if (*std::get<std::shared_ptr<const SpeciesEnum>>(const_owned) != SpeciesEnum::Water)
    throw std::runtime_error("Const rvalue variant lost its value");

  Numeric     borrowed  = 2;
  const Input read_only = borrowed;
  borrowed              = 3;
  if (*std::get<0>(read_only) != 3) throw std::runtime_error("Borrowed input copied its value");
  auto       vector       = std::make_shared<Vector>(Vector{1, 2});
  const auto const_vector = Input::from(Wsv{std::shared_ptr<Vector>(vector)});
  if (std::get<1>(const_vector).get() != vector.get()) throw std::runtime_error("Const alternative lost identity");
  auto                 pointer = std::make_shared<Numeric>(4);
  std::optional<Input> input;
  {
    Wsv owner{std::shared_ptr<Numeric>(pointer)};
    input = Input::from(owner);
    if (std::get<0>(*input).get() != pointer.get()) throw std::runtime_error("Input copied its value");
    auto output          = Output::from(owner);
    *std::get<0>(output) = 7;
    if (*pointer != 7) throw std::runtime_error("Output lost workspace identity");
    auto any = AnyInput::from(owner);
    if (std::get<std::shared_ptr<const Numeric>>(any).get() != pointer.get())
      throw std::runtime_error("Any input copied its value");
  }
  std::weak_ptr<Numeric> lifetime = pointer;
  pointer.reset();
  if (lifetime.expired() or *std::get<0>(*input) != 7) throw std::runtime_error("Input lost ownership");
  input.reset();
  if (not lifetime.expired()) throw std::runtime_error("Input leaked ownership");
  try {
    (void)Input::from(Wsv{String{"wrong type"}});
    throw std::logic_error("Unsupported alternative accepted");
  } catch (const std::runtime_error&) {}
  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return EXIT_FAILURE;
}

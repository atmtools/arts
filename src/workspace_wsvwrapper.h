#pragma once

#include <auto_wsg.h>

struct WsvValueWrapper : WsvValue {
  //! Move value into the workspace variable
  template <WorkspaceGroup T>
  WsvValueWrapper(T&& x) noexcept : WsvValue(std::make_shared<T>(std::move(x))) {}

  //! Copy value into the workspace variable
  template <WorkspaceGroup T>
  WsvValueWrapper(const T& x) : WsvValue(std::make_shared<T>(x)) {}

  //! Borrow value as workspace variable
  template <WorkspaceGroup T>
  WsvValueWrapper(T* x) noexcept : WsvValue(std::shared_ptr<T>(x, [](void*){})) {}

  //! Share value as workspace variable
  template <WorkspaceGroup T>
  WsvValueWrapper(std::shared_ptr<T>&& x) noexcept : WsvValue(std::move(x)) {}

  //! Must declare destructor to avoid incomplete type error
  WsvValueWrapper() : WsvValue(std::make_shared<Any>()) {}
  WsvValueWrapper(const WsvValueWrapper&) = default;
  WsvValueWrapper(WsvValueWrapper&&) noexcept = default;
  WsvValueWrapper& operator=(const WsvValueWrapper&) = default;
  WsvValueWrapper& operator=(WsvValueWrapper&&) noexcept = default;
};
